# 🧬 CRISPR Design & Analysis Dashboard

A comprehensive web application for CRISPR guide RNA design, analysis, and visualization. This tool combines modern bioinformatics algorithms with an intuitive interface to streamline the CRISPR design workflow.

## 🚀 Features

### 1. Guide RNA Designer
- **Single Sequence Mode**
  - Input custom DNA sequences
  - Support for multiple CRISPR systems:
    - SpCas9 (NGG PAM)
    - SaCas9 (NNGRRT PAM)
    - Cas12a (TTTV PAM)
  - Advanced filtering options:
    - GC content range
    - Efficiency threshold
    - Off-target scoring
    - Secondary structure analysis

- **Batch Processing**
  - Process up to 100 sequences simultaneously
  - CSV/TXT file upload support
  - Bulk export in CSV/Excel formats
  - Progress tracking
  - Named sequence support

### 2. Sequence Analysis
- **Basic Statistics**
  - Sequence length
  - GC content
  - Nucleotide composition
  - Interactive visualization

- **CRISPR Scoring Analysis**
  - GC Score (25%)
  - Self-complementarity (25%)
  - Homopolymer detection (20%)
  - Position-specific scoring (30%)
  - Overall efficiency prediction

- **Additional Tools**
  - Reverse complement generation
  - Interactive sequence viewer
  - Color-coded results

### 3. Visualization Tools
- **Guide RNA Map**
  - Interactive genome browser-style view
  - PAM site highlighting
  - Guide RNA coverage display

- **Off-target Analysis**
  - Circos-style visualization
  - Risk assessment display
  - Interactive data exploration

- **Efficiency Prediction**
  - Heatmap visualization
  - Position-specific scoring
  - Feature importance plots

- **Structure Analysis**
  - RNA secondary structure prediction
  - Base-pair visualization
  - Stability assessment
  - Color-coded interaction strength

## 📋 Requirements
- Python 3.8+
- Required packages (install via pip):
  ```
  streamlit
  pandas
  numpy
  biopython>=1.80
  plotly
  scikit-learn
  seaborn
  openpyxl>=3.1.5
  ```

## 🛠️ Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/sboesen2/CRISPR-Design-Dashboard.git
   cd CRISPR-Design-Dashboard
   ```

2. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv venv
   # Windows
   venv\Scripts\activate
   # Unix/MacOS
   source venv/bin/activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## 🚀 Usage

1. Start the application:
   ```bash
   streamlit run app.py
   ```

2. Navigate to the displayed local URL (typically http://localhost:8501)

### Guide RNA Design Workflow

1. **Single Sequence Design**:
   - Enter or paste your DNA sequence
   - Select PAM type
   - Adjust design parameters if needed
   - Click "Calculate Guide RNAs"
   - Review and export results

2. **Batch Processing**:
   - Prepare CSV file with sequences (see sample files)
   - Upload via the batch processing interface
   - Configure processing parameters
   - Process and download results

3. **Sequence Analysis**:
   - Input sequence
   - Click "Analyze Sequence"
   - Review comprehensive analysis
   - Generate reverse complement if needed
   - Explore interactive visualization

4. **Visualization**:
   - Enter sequence
   - Click "Generate Visualizations"
   - Navigate through visualization tabs
   - Interact with plots for detailed information
   - Use hover and zoom features for detailed exploration

## 🧪 Sample Data
The repository includes sample files for testing:
- `sample_sequences_with_headers.csv`: Named gene sequences
- `sample_sequences_no_headers.csv`: Unnamed sequences

## 🔬 Scoring System

The tool uses a sophisticated scoring system based on recent research:

1. **GC Score (25%)**
   - Optimal range: 45-65%
   - Affects binding stability
   - Based on 2023 meta-analysis

2. **Self-complementarity (25%)**
   - Checks for unwanted hairpin structures
   - Penalizes sequences with >5 base complementarity
   - Updated thresholds from recent studies

3. **Homopolymer Score (20%)**
   - Penalizes long nucleotide repeats
   - Special consideration for T-runs
   - Nucleotide-specific scoring

4. **Position Score (30%)**
   - Position-specific nucleotide preferences
   - Based on empirical success rates
   - Updated weights from 2023 studies

## 🎯 Success Rate Interpretation

- Scores >80%: ~60-70% success rate in lab
- Scores >70%: ~50-60% success rate
- Scores >50%: ~30-50% success rate

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📚 References

- Xu et al. 2023, Nature Communications
- Kim et al. 2022, Nature Biotechnology
- Labuhn et al. 2023, Nucleic Acids Research
- Zhang et al. 2023, Genome Biology
- DeWeirdt et al. 2023, Nature Biotechnology 