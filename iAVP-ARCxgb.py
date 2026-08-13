"""
Antiviral Peptide Predictor
"""
import streamlit as st
import pandas as pd
import numpy as np
from datetime import datetime
import time
import sqlite3
from pathlib import Path
import logging
import plotly.express as px
import queue
import threading
import io
import os
import sys

# Fix 1: Add current directory and parent directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

# Add multiple possible paths to sys.path
sys.path.insert(0, current_dir)  # Current directory
sys.path.insert(0, parent_dir)   # Parent directory
sys.path.insert(0, os.path.join(current_dir, 'utils'))  # utils directory
sys.path.insert(0, os.path.join(parent_dir, 'utils'))   # utils directory under parent directory

# Try to import custom functions
HAS_CUSTOM_FUNCTIONS = False
CUSTOM_MODULE_ERROR = None

try:
    # Fix 2: Try different import paths
    try:
        # Try importing from utils in current directory
        from utils.getArtFeat import get8artfeat
        st.info("✅ Import successful from utils in current directory")
    except ImportError as e1:
        # Try importing from utils in parent directory
        from utils.getArtFeat import get8artfeat
        st.info("✅ Import successful from utils in parent directory")
    
    try:
        # Try importing from ml_model in current directory
        from ml_model.iAVP_ARCfaceXGB.model_predict import deploy_predict
        st.info("✅ Import successful from ml_model in current directory")
    except ImportError as e2:
        # Try importing from ml_model in parent directory
        from ml_model.iAVP_ARCfaceXGB.model_predict import deploy_predict
        st.info("✅ Import successful from ml_model in parent directory")
    
    HAS_CUSTOM_FUNCTIONS = True
    st.success("✅ Successfully loaded custom feature extraction and model prediction functions")
    
except ImportError as e:
    HAS_CUSTOM_FUNCTIONS = False
    CUSTOM_MODULE_ERROR = str(e)
    st.warning(f"⚠️ Unable to load custom functions: {e}")
    st.warning("Will use default simplified prediction model")
    st.warning(f"Current Python path: {sys.path}")
    st.warning(f"Current working directory: {os.getcwd()}")
    st.warning(f"Script directory: {current_dir}")
    st.warning(f"Parent directory: {parent_dir}")

# Page configuration
st.set_page_config(
    page_title="Antiviral Peptide Prediction From Sequences",
    page_icon="🛡️",
    layout="wide",
    initial_sidebar_state="collapsed"  # Hide sidebar
)


# Hide Streamlit default elements
hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
section[data-testid="stSidebar"] {display: none;}
.st-emotion-cache-6q9sum {display: none;}
.st-emotion-cache-1wbqy5l {display: none;}

/* Disable button click states */
button:disabled {
    opacity: 0.6;
    cursor: not-allowed;
}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/antiviral_peptide.log'),
    ]
)
logger = logging.getLogger(__name__)

# Initialize Session State variables
if 'prediction_in_progress' not in st.session_state:
    st.session_state.prediction_in_progress = False
if 'prediction_results' not in st.session_state:
    st.session_state.prediction_results = None
if 'prediction_errors' not in st.session_state:
    st.session_state.prediction_errors = []
if 'prediction_progress' not in st.session_state:
    st.session_state.prediction_progress = 0.0
if 'prediction_status' not in st.session_state:
    st.session_state.prediction_status = ""
if 'prediction_substatus' not in st.session_state:
    st.session_state.prediction_substatus = ""
if 'prediction_thread' not in st.session_state:
    st.session_state.prediction_thread = None
if 'progress_queue' not in st.session_state:
    st.session_state.progress_queue = None
if 'sequences_to_process' not in st.session_state:
    st.session_state.sequences_to_process = []
if 'reset_requested' not in st.session_state:
    st.session_state.reset_requested = False
if 'fasta_input' not in st.session_state:
    st.session_state.fasta_input = ""
if 'uploaded_file_key' not in st.session_state:
    st.session_state.uploaded_file_key = 0
if 'background_results' not in st.session_state:
    st.session_state.background_results = None
if 'background_errors' not in st.session_state:
    st.session_state.background_errors = []
if 'background_finished' not in st.session_state:
    st.session_state.background_finished = False
if 'batch_processing' not in st.session_state:
    st.session_state.batch_processing = False
if 'current_batch' not in st.session_state:
    st.session_state.current_batch = 0
if 'total_batches' not in st.session_state:
    st.session_state.total_batches = 1
if 'batch_results' not in st.session_state:
    st.session_state.batch_results = []
if 'sequence_limit_exceeded' not in st.session_state:
    st.session_state.sequence_limit_exceeded = False
if 'original_sequence_count' not in st.session_state:
    st.session_state.original_sequence_count = 0
if 'seq_pd' not in st.session_state:
    st.session_state.seq_pd = None
if 'feature_pd' not in st.session_state:
    st.session_state.feature_pd = None
if 'all_features_extracted' not in st.session_state:
    st.session_state.all_features_extracted = False
if 'current_stage' not in st.session_state:
    st.session_state.current_stage = ""
if 'stage_progress' not in st.session_state:
    st.session_state.stage_progress = 0.0

# Database access logging
def log_visit(ip, page):
    try:
        Path("data").mkdir(exist_ok=True)
        conn = sqlite3.connect("data/visits.db", check_same_thread=False)
        cursor = conn.cursor()
        cursor.execute(
            'CREATE TABLE IF NOT EXISTS visits (id INTEGER PRIMARY KEY AUTOINCREMENT, ip TEXT, page TEXT, timestamp DATETIME DEFAULT CURRENT_TIMESTAMP)'
        )
        cursor.execute(
            'INSERT INTO visits (ip, page) VALUES (?, ?)',
            (ip, page)
        )
        conn.commit()
    except Exception as e:
        logger.error(f"Failed to log visit: {e}")

def get_client_ip():
    """Get client IP address"""
    try:
        if hasattr(st, 'context') and hasattr(st.context, 'headers'):
            headers = st.context.headers
            if 'X-Forwarded-For' in headers:
                return headers['X-Forwarded-For'].split(',')[0].strip()
    except:
        pass
    return "127.0.0.1"

# Log visit
client_ip = get_client_ip()
log_visit(client_ip, "antiviral_peptide")

# Page title
st.markdown("""
<div style="text-align: center; margin-bottom: 2rem;">
    <h1 style="color: #1a73e8; font-size: 2.5rem; margin-bottom: 0.5rem;">🛡️ Antiviral Peptide Predictor</h1>
    <p style="color: #5f6368; font-size: 1.1rem;">Predict antiviral activity of peptide sequences</p>
</div>
""", unsafe_allow_html=True)

# Peptide sequence processing functions
def parse_fasta_to_dataframe(text):
    """
    Parse FASTA format text and convert to DataFrame
    
    Parameters:
    text: FASTA format text
    
    Returns:
    pd.DataFrame: DataFrame containing ID and Sequence columns
    """
    sequences = []
    current_id = ""
    current_seq = ""
    
    for line in text.split('\n'):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('>'):
            if current_id and current_seq:
                sequences.append({'ID': current_id, 'Sequence': current_seq})
            current_id = line[1:].split()[0] if ' ' in line[1:] else line[1:]
            current_seq = ""
        else:
            current_seq += line.upper()
    
    if current_id and current_seq:
        sequences.append({'ID': current_id, 'Sequence': current_seq})
    
    # Convert to DataFrame
    seq_pd = pd.DataFrame(sequences)
    return seq_pd

def validate_peptide(sequence):
    """Validate peptide sequence"""
    valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
    seq_upper = sequence.upper()
    
    if not sequence:
        return False, "Sequence is empty"
    
    if len(sequence) < 3:
        return False, "Sequence too short (min 3 amino acids)"
    
    if len(sequence) > 100:
        return False, "Sequence too long (max 100 amino acids)"
    
    invalid = set(seq_upper) - valid_aas
    if invalid:
        return False, f"Invalid amino acids: {', '.join(invalid)}"
    
    return True, "Valid"

def predict_antiviral_simple(sequence):
    """Simplified antiviral peptide prediction (fallback solution)"""
    seq = sequence.upper()
    
    # Calculate features
    cys_content = seq.count('C') / len(seq)
    aromatic_content = (seq.count('F') + seq.count('W') + seq.count('Y')) / len(seq)
    basic_content = (seq.count('K') + seq.count('R') + seq.count('H')) / len(seq)
    
    # Base score
    score = cys_content * 0.4 + aromatic_content * 0.3 + basic_content * 0.3
    
    # Length adjustment
    length = len(seq)
    if 5 <= length <= 30:
        length_factor = 1.0
    else:
        length_factor = 0.7
    
    final_score = min(1.0, score * length_factor)
    
    # Determine prediction
    if final_score > 0.6:
        label = "Antiviral"
        probability = final_score
    else:
        label = "Non-antiviral"
        probability = 1.0 - final_score
    
    # Confidence
    if probability > 0.8:
        confidence = "High"
    elif probability > 0.6:
        confidence = "Medium"
    else:
        confidence = "Low"
    
    return {
        'label': label,
        'probability': probability,
        'confidence': confidence,
        'score': final_score
    }

def process_sequences_background(seq_pd, progress_callback=None):
    """
    Background sequence processing - Modified to: extract features one by one, then predict one by one after all features are extracted
    
    Parameters:
    seq_pd: DataFrame, containing ID and Sequence columns
    progress_callback: Progress callback function
    
    Returns:
    tuple: (results, errors)
    """
    results = []
    errors = []
    all_features_list = []  # Store all extracted features
    
    # Stage 1: Check sequence count
    if progress_callback:
        progress_callback(0.0, "Starting sequence processing...", "Initialization")
    
    # Check sequence count, if exceeds 200, only process first 200
    if len(seq_pd) > 200:
        processed_seq_pd = seq_pd.head(200)
        
        # Set flag
        st.session_state.sequence_limit_exceeded = True
        st.session_state.original_sequence_count = len(seq_pd)
        
        if progress_callback:
            progress_callback(0.05, f"Detected {len(seq_pd)} sequences, exceeding 200 limit, will only process first 200", "Initialization")
    else:
        processed_seq_pd = seq_pd
        st.session_state.sequence_limit_exceeded = False
        st.session_state.original_sequence_count = len(seq_pd)
    
    total_sequences = len(processed_seq_pd)
    
    if total_sequences == 0:
        return results, errors
    
    # Stage 2: Extract features one by one
    if progress_callback:
        progress_callback(0.1, f"Starting feature extraction one by one, total {total_sequences} sequences", "Feature Extraction")
    
    for i, (_, row) in enumerate(processed_seq_pd.iterrows()):
        seq_id = row['ID']
        sequence = row['Sequence']
        
        # Update feature extraction progress
        feature_progress = 0.1 + (i / total_sequences) * 0.4
        if progress_callback:
            progress_callback(feature_progress, f"Feature extraction progress: {i+1}/{total_sequences} - {seq_id}", "Feature Extraction")
        
        try:
            if HAS_CUSTOM_FUNCTIONS:
                # Extract features for single sequence
                single_seq_pd = pd.DataFrame([{'ID': seq_id, 'Sequence': sequence}])
                features = get8artfeat(single_seq_pd)
                
                # Ensure features is DataFrame format
                if not isinstance(features, pd.DataFrame):
                    raise ValueError(f"Feature extraction function returned wrong type: {type(features)}")
                
                # Save feature data
                all_features_list.append({
                    'ID': seq_id,
                    'Sequence': sequence,
                    'features': features
                })
            else:
                # Simplified model doesn't need feature extraction
                all_features_list.append({
                    'ID': seq_id,
                    'Sequence': sequence,
                    'features': None
                })
                
        except Exception as e:
            error_msg = f"Sequence {seq_id} feature extraction failed: {str(e)}"
            errors.append(error_msg)
            logger.error(error_msg)
            continue
        
        # Simulate processing time
        time.sleep(0.05)
    
    if progress_callback:
        progress_callback(0.5, f"Feature extraction completed, extracted features for {len(all_features_list)} sequences", "Feature Extraction Complete")
    
    # Stage 3: Machine learning prediction one by one
    if progress_callback:
        progress_callback(0.5, f"Starting machine learning prediction one by one, total {len(all_features_list)} sequences", "Model Prediction")
    
    for i, feature_item in enumerate(all_features_list):
        seq_id = feature_item['ID']
        sequence = feature_item['Sequence']
        features = feature_item['features']
        
        # Update prediction progress
        predict_progress = 0.5 + (i / len(all_features_list)) * 0.4
        if progress_callback:
            progress_callback(predict_progress, f"Model prediction progress: {i+1}/{len(all_features_list)} - {seq_id}", "Model Prediction")
        
        try:
            if HAS_CUSTOM_FUNCTIONS and features is not None:
                # Predict one by one
                prediction = deploy_predict(features)
                
                # Ensure prediction result is DataFrame format
                if not isinstance(prediction, pd.DataFrame):
                    if isinstance(prediction, dict):
                        prediction = pd.DataFrame([prediction])
                    elif isinstance(prediction, list):
                        prediction = pd.DataFrame(prediction)
                    else:
                        raise ValueError("Model prediction function returned wrong format")
                
                if len(prediction) > 0:
                    pred_row = prediction.iloc[0]
                    label = pred_row.get('label', 'Unknown')
                    probability = pred_row.get('probability', 0.5)
                else:
                    raise ValueError("Prediction result is empty")
                    
            else:
                # Use simplified model
                simple_pred = predict_antiviral_simple(sequence)
                label = simple_pred['label']
                probability = simple_pred['probability']
                confidence = simple_pred['confidence']
                score = simple_pred['score']
            
            # Save result
            results.append({
                'ID': seq_id,
                'Sequence': sequence,
                'Length': len(sequence),
                'Prediction': label,
                'Probability(%)': f"{probability*100:.3f}",
            })
            
        except Exception as e:
            error_msg = f"Sequence {seq_id} prediction failed: {str(e)}"
            errors.append(error_msg)
            logger.error(error_msg)
            continue
        
        # Simulate processing time
        time.sleep(0.05)
    
    # Stage 4: Complete processing
    if progress_callback:
        progress_callback(0.9, "Organizing and formatting results...", "Result Organization")
        time.sleep(0.5)
        progress_callback(1.0, "Processing completed!", "Complete")
    
    return results, errors

def start_prediction(seq_pd):
    """Start prediction - using thread-safe queue"""
    st.session_state.prediction_in_progress = True
    st.session_state.prediction_results = None
    st.session_state.prediction_errors = []
    st.session_state.prediction_progress = 0.0
    st.session_state.prediction_status = "Starting prediction..."
    st.session_state.prediction_substatus = ""
    st.session_state.current_stage = ""
    st.session_state.stage_progress = 0.0
    st.session_state.background_finished = False
    st.session_state.background_results = None
    st.session_state.background_errors = []
    st.session_state.batch_results = []
    st.session_state.seq_pd = seq_pd
    st.session_state.feature_pd = None
    st.session_state.all_features_extracted = False
    
    # No longer using batch processing since we're processing one by one
    st.session_state.batch_processing = False
    st.session_state.current_batch = 0
    st.session_state.total_batches = 1
    
    # Create queue for thread communication
    st.session_state.progress_queue = queue.Queue()
    
    def background_worker(seq_df, q):
        """Background worker thread"""
        def progress_callback(progress, status, stage=""):
            q.put(('progress', progress, status, stage))
        
        try:
            results, errors = process_sequences_background(seq_df, progress_callback)
            q.put(('results', results, errors))
        except Exception as e:
            q.put(('error', str(e)))
    
    # Start background thread
    thread = threading.Thread(
        target=background_worker,
        args=(seq_pd, st.session_state.progress_queue),
        daemon=True
    )
    thread.start()
    st.session_state.prediction_thread = thread
    
    # Single click execution, no need for additional clicks
    st.rerun()

def reset_application():
    """Reset application"""
    st.session_state.reset_requested = True
    st.rerun()

# Handle reset request
if st.session_state.reset_requested:
    st.session_state.prediction_in_progress = False
    st.session_state.prediction_results = None
    st.session_state.prediction_errors = []
    st.session_state.prediction_progress = 0.0
    st.session_state.prediction_status = ""
    st.session_state.prediction_substatus = ""
    st.session_state.current_stage = ""
    st.session_state.stage_progress = 0.0
    st.session_state.prediction_thread = None
    st.session_state.progress_queue = None
    st.session_state.sequences_to_process = []
    st.session_state.fasta_input = ""
    st.session_state.seq_pd = None
    st.session_state.feature_pd = None
    st.session_state.all_features_extracted = False
    st.session_state.background_results = None
    st.session_state.background_errors = []
    st.session_state.background_finished = False
    st.session_state.batch_processing = False
    st.session_state.current_batch = 0
    st.session_state.total_batches = 1
    st.session_state.batch_results = []
    st.session_state.sequence_limit_exceeded = False
    st.session_state.original_sequence_count = 0
    st.session_state.uploaded_file_key += 1
    st.session_state.reset_requested = False
    st.rerun()

# Handle progress updates
if st.session_state.prediction_in_progress and st.session_state.progress_queue is not None:
    try:
        # Process all messages in queue
        while not st.session_state.progress_queue.empty():
            try:
                msg = st.session_state.progress_queue.get_nowait()
                
                if msg[0] == 'progress':
                    _, progress, status, stage = msg
                    st.session_state.prediction_progress = progress
                    st.session_state.prediction_status = status
                    if stage:
                        st.session_state.current_stage = stage
                elif msg[0] == 'results':
                    _, results, errors = msg
                    st.session_state.background_results = results
                    st.session_state.background_errors = errors
                    st.session_state.background_finished = True
                elif msg[0] == 'error':
                    _, error_msg = msg
                    st.session_state.prediction_errors = [f"Background processing error: {error_msg}"]
                    st.session_state.background_finished = True
            except queue.Empty:
                break
    except Exception as e:
        st.session_state.prediction_errors = [f"Error processing queue: {str(e)}"]
        st.session_state.background_finished = True

# Check if background processing is complete
if st.session_state.background_finished and st.session_state.prediction_in_progress:
    st.session_state.prediction_in_progress = False
    
    # Process results
    if st.session_state.background_results:
        results_df = pd.DataFrame(st.session_state.background_results)
        st.session_state.prediction_results = results_df
    else:
        st.session_state.prediction_results = pd.DataFrame()
    
    st.session_state.prediction_errors = st.session_state.background_errors
    
    # Reset background state
    st.session_state.background_results = None
    st.session_state.background_errors = []
    st.session_state.background_finished = False
    
    st.rerun()

# Display progress
if st.session_state.prediction_in_progress:
    st.markdown("### ⏳ Processing Status")
    
    # Display current stage
    if st.session_state.current_stage:
        st.info(f"**Current Stage:** {st.session_state.current_stage}")
    
    # Main progress bar
    progress_bar = st.progress(st.session_state.prediction_progress)
    
    # Status text
    if st.session_state.prediction_status:
        st.info(st.session_state.prediction_status)
    
    # Display detailed progress stages
    with st.expander("🔍 Detailed Progress", expanded=False):
        stages = [
            ("Initialization", 0.0, 0.1),
            ("Feature Extraction", 0.1, 0.5),
            ("Model Prediction", 0.5, 0.9),
            ("Result Organization", 0.9, 1.0)
        ]
        
        for stage_name, start_progress, end_progress in stages:
            if st.session_state.prediction_progress >= start_progress:
                stage_progress = min(1.0, (st.session_state.prediction_progress - start_progress) / (end_progress - start_progress))
                st.progress(stage_progress, text=f"{stage_name}: {stage_progress*100:.1f}%")
    
    # Add polling
    time.sleep(0.1)
    st.rerun()

# Input section
st.markdown("### 📥 Input Peptide Sequences")

input_method = st.radio(
    "Select input method:",
    ["Text Input", "File Upload"],
    horizontal=True,
    key="input_method_antiviral"
)

sequences = []
input_content = ""
seq_pd = None

if input_method == "Text Input":
    if 'fasta_input' in st.session_state:
        input_content = st.session_state.fasta_input
    
    input_content = st.text_area(
        "Enter peptide sequences:",
        value=input_content,
        height=200,
        placeholder=">peptide1\nACDEFGHIK\n>peptide2\nCCHHMMWWYY",
        key=f"fasta_input_area_{st.session_state.uploaded_file_key}",
        disabled=st.session_state.prediction_in_progress
    )
    
    st.session_state.fasta_input = input_content
    
    if input_content and not st.session_state.prediction_in_progress:
        # Parse to DataFrame
        seq_pd = parse_fasta_to_dataframe(input_content)
        
        if not seq_pd.empty:
            st.info(f"✅ Found {len(seq_pd)} sequences")
            
            # Display sequence limit warning
            if len(seq_pd) > 200:
                st.warning(f"⚠️ Detected {len(seq_pd)} sequences, exceeding 200 limit, will only process first 200.")
            
            # Display sequence preview
            with st.expander("🔍 Preview Sequences"):
                st.dataframe(seq_pd.head(10),width='stretch')
                if len(seq_pd) > 10:
                    st.write(f"Showing first 10 sequences, total {len(seq_pd)} sequences")
else:
    uploaded_file = st.file_uploader(
        "Upload FASTA file:",
        type=['fasta', 'fa', 'txt', 'pep'],
        key=f"file_upload_antiviral_{st.session_state.uploaded_file_key}",
        disabled=st.session_state.prediction_in_progress
    )
    
    if uploaded_file is not None and not st.session_state.prediction_in_progress:
        try:
            input_content = uploaded_file.getvalue().decode('utf-8')
            # Parse to DataFrame
            seq_pd = parse_fasta_to_dataframe(input_content)
            st.session_state.fasta_input = input_content
            
            if not seq_pd.empty:
                st.success(f"✅ Successfully loaded {len(seq_pd)} sequences")
                
                # Display sequence limit warning
                if len(seq_pd) > 200:
                    st.warning(f"⚠️ Detected {len(seq_pd)} sequences, exceeding 200 limit, will only process first 200.")
                
                # Display sequence preview
                with st.expander("🔍 Preview Sequences"):
                    st.dataframe(seq_pd.head(10),width='stretch')
                    if len(seq_pd) > 10:
                        st.write(f"Showing first 10 sequences, total {len(seq_pd)} sequences")
        except Exception as e:
            st.error(f"Error reading file: {e}")

# Submit button and processing
if seq_pd is not None and not seq_pd.empty and not st.session_state.prediction_in_progress:
    col1, col2 = st.columns([3, 1])
    
    with col1:
        # Single click to execute prediction
        if st.button("🧠 Predict Antiviral Activity", type="primary",width='stretch'):
            start_prediction(seq_pd)
    
    with col2:
        # Single click to reset
        if st.button("🔄 Reset", type="secondary",width='stretch'):
            reset_application()

# Display results
if st.session_state.prediction_results is not None and not st.session_state.prediction_in_progress:
    results_df = st.session_state.prediction_results
    errors = st.session_state.prediction_errors
    
    # Display sequence limit notification
    if st.session_state.sequence_limit_exceeded:
        st.warning(f"""
        ⚠️ **Sequence Limit Notification**
        
        You entered {st.session_state.original_sequence_count} sequences, exceeding the 200 sequence limit.
        The system has automatically processed the first 200 sequences, the remaining sequences were not processed.
        
        To process all sequences, please submit in batches.
        """)
    
    # Display errors
    if errors:
        st.warning("⚠️ Some sequences had issues:")
        for error in errors:
            st.error(error)
    
    # Display prediction results
    if not results_df.empty:
        st.markdown("### 📊 Results")
        
        st.dataframe(
            results_df,
            column_config={
                "ID": "Sequence ID",
                "Sequence": "Amino Acid Sequence",
                "Length": "Length",
                "Prediction": "Prediction",
                "Probability": "Probability"
            },
            hide_index=True,
           width='stretch'
        )
        
        # Visualization only (removed statistical metrics as requested)
        st.markdown("### 📈 Visualization")
        
        col_viz1, col_viz2 = st.columns(2)
        
        with col_viz1:
            # Prediction distribution pie chart
            fig1 = px.pie(
                results_df,
                names='Prediction',
                title='Prediction Distribution',
                color='Prediction',
                color_discrete_map={
                    'Antiviral': '#4caf50',
                    'Non-antiviral': '#f44336'
                }
            )
            fig1.update_layout(plot_bgcolor='white', paper_bgcolor='white')
            st.plotly_chart(fig1,width='stretch')
        
        with col_viz2:
            # Probability histogram
            results_df['Probability_num'] = results_df['Probability(%)'].astype(float)
            fig2 = px.histogram(
                results_df,
                x='Probability_num',
                color='Prediction',
                title='Probability Distribution',
                nbins=20
            )
            fig2.update_layout(plot_bgcolor='white', paper_bgcolor='white')
            st.plotly_chart(fig2,width='stretch')
        
        # Export options
        st.markdown("### 💾 Export Results")
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        col_export1, col_export2, col_export3 = st.columns(3)
        
        with col_export1:
            csv_data = results_df.drop(['Probability_num'], axis=1, errors='ignore').to_csv(index=False)
            st.download_button(
                label="📥 Download CSV",
                data=csv_data,
                file_name=f"antiviral_predictions_{timestamp}.csv",
                mime="text/csv",
               width='stretch'
            )
        
        with col_export2:
            json_data = results_df.drop(['Probability_num'], axis=1, errors='ignore').to_json(orient='records', indent=2)
            st.download_button(
                label="📥 Download JSON",
                data=json_data,
                file_name=f"antiviral_predictions_{timestamp}.json",
                mime="application/json",
               width='stretch'
            )
        
        with col_export3:
            # Create Excel file
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                # Prediction results
                results_df.drop(['Probability_num'], axis=1, errors='ignore').to_excel(
                    writer, index=False, sheet_name='Predictions'
                )
                
                # Add summary information
                summary_data = {
                    'Metric': ['Total Sequences', 'Antiviral', 'Non-antiviral', 'Average Probability', 'Processing Date'],
                    'Value': [
                        len(results_df),
                        len(results_df[results_df['Prediction'] == 'Antiviral']),
                        len(results_df[results_df['Prediction'] == 'Non-antiviral']),
                        results_df['Probability(%)'].astype(float).mean(),
                        datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    ]
                }
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, index=False, sheet_name='Summary')
            
            st.download_button(
                label="📥 Download Excel",
                data=output.getvalue(),
                file_name=f"antiviral_predictions_{timestamp}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
               width='stretch'
            )
        
        # New prediction and reset buttons
        st.markdown("---")
        col_new1, col_new2 = st.columns(2)
        with col_new1:
            # Single click for new prediction
            if st.button("🧠 Run New Prediction", type="primary",width='stretch'):
                st.session_state.prediction_results = None
                st.session_state.prediction_errors = []
                st.session_state.seq_pd = None
                st.session_state.feature_pd = None
                st.session_state.sequence_limit_exceeded = False
                st.rerun()
        with col_new2:
            # Single click to reset all
            if st.button("🔄 Reset All", type="secondary",width='stretch'):
                reset_application()

# Model information
with st.expander("ℹ️ Model Information"):
    if HAS_CUSTOM_FUNCTIONS:
        st.markdown("""
        ### 🧬 Machine Learning Model
        
        This tool uses a trained machine learning model for antiviral peptide prediction:
        
        **Workflow:**
        1. **FASTA Parsing**: Parse input sequences to DataFrame (seq_pd)
        2. **Feature Extraction**: Extract features using `get8artfeat(seq_pd)` -> feature_pd
        3. **Model Prediction**: Predict with `deploy_predict(feature_pd)`
        4. **Result Display**: Show predictions with probabilities
        
        **Data Flow:**
        - Input: FASTA format -> seq_pd (DataFrame with ID, Sequence)
        - Feature Extraction: seq_pd -> feature_pd (DataFrame with ID, Sequence, features)
        - Prediction: feature_pd -> prediction_results (DataFrame)
        
        **Processing Strategy:**
        - **Sequential Feature Extraction**: Extract features for all sequences first
        - **Sequential Prediction**: Make predictions for all sequences after feature extraction
        - Maximum 200 sequences per submission
        
        **Progress Stages:**
        1. Initialization
        2. Feature Extraction
        3. Model Prediction
        4. Results Consolidation
        
        **Note:** Ensure your utils/ and ml_model/ directories contain the required modules.
        """)
    else:
        st.markdown("""
        ### ⚠️ Simplified Model (Fallback)
        
        **Note:** Custom ML modules not found. Using simplified rule-based model.
        
        To use the ML model, add these files:
        1. `utils/getArtFeat.py` - Feature extraction function
        2. `ml_model/iAVP_ARCfaceXGB/model_predict.py` - Model prediction function (deploy_predict)
        
        **Current simplified model uses:**
        1. Cysteine content (40% weight)
        2. Aromatic amino acids (30% weight)
        3. Basic amino acids (30% weight)
        4. Sequence length adjustment
        
        **Processing Strategy:**
        - **Sequential Feature Extraction**: Extract features for all sequences first
        - **Sequential Prediction**: Make predictions for all sequences after feature extraction
        - Maximum 200 sequences per submission
        """)

# Footer button area
st.markdown("---")
footer_col1, footer_col2 = st.columns(2)

with footer_col1:
    # Single click to return to main page
    if st.button("← Back to Main", type="secondary",width='stretch'):
        try:
            st.switch_page("app.py")
        except Exception as e:
            st.rerun()

with footer_col2:
    # Single click to reset application
    if st.button("🔄 Reset Application", type="secondary",width='stretch'):
        reset_application()

# Add ICP备案信息在页面底部
st.markdown("""
<div style="text-align: center; margin-top: 3rem; padding: 1rem; border-top: 1px solid #e0e0e0;">
    <p style="color: #666; font-size: 0.9rem;">
        <a href="https://beian.miit.gov.cn/" target="_blank" style="text-decoration: none; color: #666;">
            蜀ICP备19038779号-1
        </a>
    </p>
</div>
""", unsafe_allow_html=True)

# Standalone run check
if __name__ == "__main__":
    st.info("✅ Application is running in standalone mode")