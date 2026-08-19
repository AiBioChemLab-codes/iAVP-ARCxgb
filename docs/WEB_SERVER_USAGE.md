# Web-server usage

## Online server

Open:

https://iservers.aibiochem.net/iAVP-ARCxgb

## Input format

Paste FASTA-formatted peptide sequences or upload a FASTA/text file.

```text
>peptide_1
IFWDCWAPEEPACQDFLGAMIH
>peptide_2
LSELDDRADALQAGASQFETSAAKLKRKYWWKN
```

Input requirements:

- Standard amino-acid characters only: `ACDEFGHIKLMNPQRSTVWY`
- Length: 3-50 amino acids
- Up to 200 sequences per submission

## Output

For each sequence, the server reports:

- sequence ID
- peptide sequence
- sequence length
- predicted class
- AVP-class probability

Results can be exported as CSV, JSON, or Excel files.

## Local deployment

```bash
pip install -r requirements.txt
streamlit run iAVP-ARCxgb.py
```

The Streamlit application uses the serialized deployment model under `ml_model/iAVP_ARCfaceXGB/`.
