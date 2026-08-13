import joblib
import torch
import os
 
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import sys
 


import warnings
warnings.filterwarnings('ignore')

# ===================== GPU显存自动管理 =====================
def clear_gpu_memory():
    """清理GPU显存"""
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()

class TabularEmbedNet(nn.Module):
    def __init__(self, input_dim, embed_dim=64):
        super().__init__()
        self.input_dim = input_dim
        self.embed_dim = embed_dim
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256), nn.BatchNorm1d(256), nn.LeakyReLU(0.2), nn.Dropout(0.3),
            nn.Linear(256, 128), nn.BatchNorm1d(128), nn.LeakyReLU(0.2), nn.Dropout(0.3),
            nn.Linear(128, embed_dim)
        )
    
    def forward(self, x):
        return self.net(x)
    
    def save_model(self, path):
        """保存完整模型（包括结构和参数）"""
        torch.save({
            'model_state_dict': self.state_dict(),
            'input_dim': self.input_dim,
            'embed_dim': self.embed_dim
        }, path)
    
    @classmethod
    def load_model(cls, path, device=None):
        """加载完整模型"""
        if device is None:
            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        checkpoint = torch.load(path, map_location=device)
        model = cls(checkpoint['input_dim'], checkpoint['embed_dim'])
        model.load_state_dict(checkpoint['model_state_dict'])
        model.to(device)
        model.eval()
        return model
    

def deploy_predict(feature_pd, gpu_id=-1):
    

    xgb_deploy_path=os.path.join(os.path.dirname(__file__),"iAVP_arcfaceXGB_deploy_train_init_xgb_deploy.joblib")
    arcface_deploy_path=os.path.join(os.path.dirname(__file__),"iAVP_arcfaceXGB_deploy_train_init_arcface_deploy.pth")
    scaler_deploy_path=os.path.join(os.path.dirname(__file__),"iAVP_arcfaceXGB_deploy_train_init_scaler_deploy.joblib")
    params_deploy_path=os.path.join(os.path.dirname(__file__),"iAVP_arcfaceXGB_deploy_train_init_best_params.csv")
    scaler_init_path=os.path.join(os.path.dirname(__file__),"iAVP_arcfaceXGB_deploy_train_init_scaler.joblib")

    
    deploy_paths={

        "xgb_path": xgb_deploy_path,
        "arcface_path": arcface_deploy_path,
        "scaler_path": scaler_deploy_path,
        "params_path": params_deploy_path,
        "scaler_init_path": scaler_init_path


    }
    

    
    # Setup Device
    if torch.cuda.is_available() and gpu_id >= 0:
        device = torch.device(f'cuda:{gpu_id}')
    else:
        device = torch.device('cpu')
    
    print(f"[部署预测] 使用设备: {device}")
    
    # 1. Load Model
    print(f"[部署预测] 加载部署模型...")
    xgb_model = joblib.load(deploy_paths["xgb_path"])
    arc_model = TabularEmbedNet.load_model(deploy_paths["arcface_path"], device=device)
     
    scaler = joblib.load(deploy_paths["scaler_path"])
    params_df = pd.read_csv(deploy_paths["params_path"])
    best_th = params_df.iloc[0]["optimal_threshold"]
    
    # 2. 加载测试数据
    df_test = feature_pd
    id_test = df_test.iloc[:, 0].values
    id_test_seqs = df_test.iloc[:, 1].values
    X_test_original = df_test.iloc[:, 2:].values

    id_test = df_test.iloc[:, 0].values
   
    X_test_original = df_test.iloc[:, 2:].values
    
    # 3. 标准化
    X_test_scaled = scaler.transform(X_test_original)
    
    # 4. 生成嵌入特征
    X_tensor = torch.tensor(X_test_scaled, dtype=torch.float32).to(device)
    with torch.no_grad():
        arc_model.eval()
        emb_test = arc_model(X_tensor).cpu().numpy()
    
    # 5. 合并特征
    X_test_combined = np.hstack([X_test_scaled, emb_test])
    
    # 6. 预测
    y_prob = xgb_model.predict_proba(X_test_combined)[:, 1]
    y_pred = (y_prob >= best_th).astype(int)
    
    # 7. 计算指标
   
    # 生成时间戳
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
   # 8.返回预测结果
    pred_df = pd.DataFrame({
        "ID": id_test,
        "Sequence": id_test_seqs,
        "label": y_pred,
        "probability": y_prob
    })
     
    #pred_df.to_csv("pred.csv")
    
    # 清理GPU内存
    clear_gpu_memory()
    
    return pred_df

