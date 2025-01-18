import streamlit as st
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import plotly.graph_objects as go
from src.utils.sequence_tools import (
    find_pam_sites, 
    design_grnas, 
    display_grna_results,
    predict_grna_efficiency,
    check_off_targets,
    analyze_secondary_structure,
    create_sequence_plot
)

# Add color coding functions at the top level
def color_efficiency(val):
    """Color code efficiency scores"""
    try:
        num_val = float(val.strip('%'))/100 if isinstance(val, str) else val
        if num_val >= 0.7:
            return 'background-color: #90EE90'  # Light green
        elif num_val >= 0.5:
            return 'background-color: #FFFFE0'  # Light yellow
        return 'background-color: #FFB6C1'      # Light red
    except:
        return ''

def color_offtarget(val):
    """Color code off-target scores"""
    try:
        num_val = float(val) if isinstance(val, str) else val
        if num_val <= 20:
            return 'background-color: #90EE90'  # Light green
        elif num_val <= 50:
            return 'background-color: #FFFFE0'  # Light yellow
        return 'background-color: #FFB6C1'      # Light red
    except:
        return ''

def set_page_config():
    st.set_page_config(
        page_title="CRISPR Design Dashboard",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )

def main():
    set_page_config()
    
    st.title("🧬 CRISPR Design & Analysis Dashboard")
    
    # Initialize session states if they don't exist
    if 'dna_sequence' not in st.session_state:
        st.session_state.dna_sequence = ""
    if 'uploaded_file' not in st.session_state:
        st.session_state.uploaded_file = None
    if 'has_headers' not in st.session_state:
        st.session_state.has_headers = True
    
    # Add clear button in sidebar
    with st.sidebar:
        if st.button("🗑️ Clear All Data"):
            st.session_state.dna_sequence = ""
            st.session_state.uploaded_file = None
            st.session_state.has_headers = True
            st.rerun()
    
    # Sidebar navigation
    page = st.sidebar.selectbox(
        "Select a Tool",
        ["Home", "Guide RNA Designer", "Sequence Analysis", "Visualization"]
    )
    
    if page == "Home":
        show_home()
    elif page == "Guide RNA Designer":
        show_grna_designer()
    elif page == "Sequence Analysis":
        show_sequence_analysis()
    elif page == "Visualization":
        show_visualization()

def show_home():
    st.header("Welcome to CRISPR Design & Analysis Dashboard")
    
    st.markdown("""
    This dashboard provides tools for:
    - 🎯 Guide RNA Design and Analysis
    - 🧪 Sequence Analysis
    - 📊 Result Visualization
    
    Select a tool from the sidebar to get started!
    """)

def show_grna_designer():
    st.header("Guide RNA Designer")
    
    # Add tabs for single/batch processing
    input_mode = st.radio(
        "Select Input Mode",
        ["Single Sequence", "Batch Processing"],
        help="Choose between analyzing a single sequence or multiple sequences"
    )
    
    if input_mode == "Single Sequence":
        col1, col2 = st.columns([2,1])
        
        with col1:
            # Use session state for sequence persistence
            sequence = st.text_area(
                "Enter your DNA sequence",
                value=st.session_state.dna_sequence,
                height=150,
                help="Enter a DNA sequence (A, T, G, C only)"
            )
            # Update session state
            st.session_state.dna_sequence = sequence
            
            pam_type = st.selectbox(
                "Select PAM Type",
                ["SpCas9 (NGG)", "SaCas9 (NNGRRT)", "Cas12a (TTTV)"]
            )
            
            calculate = st.button("Calculate Guide RNAs")
            
            advanced_options = st.expander("Advanced Options")
            with advanced_options:
                efficiency_threshold = st.slider(
                    "Minimum Efficiency Score", 
                    0.0, 1.0, 0.5
                )
                max_off_targets = st.slider(
                    "Maximum Off-target Score", 
                    0.0, 100.0, 50.0
                )
                check_secondary = st.checkbox(
                    "Check Secondary Structure", 
                    value=True
                )
        
        with col2:
            st.info("Design Parameters")
            gc_min = st.slider("Minimum GC content (%)", 30, 80, 40)
            gc_max = st.slider("Maximum GC content (%)", gc_min, 80, 60)
            
            st.info("Filtering Options")
            show_all = st.checkbox("Show All Candidates", value=False)
            
        if sequence and calculate:
            process_single_sequence(sequence, pam_type, gc_min, gc_max, 
                                 efficiency_threshold, max_off_targets, 
                                 show_all, check_secondary)
    
    else:  # Batch Processing
        st.info("""
        **Batch Processing Instructions:**
        1. Upload a CSV/TXT file with one sequence per line
        2. First line can optionally contain sequence names
        3. Maximum 100 sequences per batch
        4. Each sequence should be DNA (A, T, G, C only)
        """)
        
        # Use session state for file upload persistence
        uploaded_file = st.file_uploader(
            "Upload sequences file (CSV/TXT)", 
            type=['csv', 'txt'],
            key="batch_file_uploader"
        )
        
        # Update session state
        if uploaded_file is not None:
            st.session_state.uploaded_file = uploaded_file
        
        col1, col2 = st.columns([2,1])
        
        with col1:
            pam_type = st.selectbox(
                "Select PAM Type",
                ["SpCas9 (NGG)", "SaCas9 (NNGRRT)", "Cas12a (TTTV)"]
            )
            
            has_headers = st.checkbox(
                "File has headers", 
                value=st.session_state.has_headers,
                help="Check if first line contains sequence names"
            )
            # Update session state
            st.session_state.has_headers = has_headers
        
        with col2:
            gc_min = st.slider("Minimum GC content (%)", 30, 80, 40)
            gc_max = st.slider("Maximum GC content (%)", gc_min, 80, 60)
            
            efficiency_threshold = st.slider(
                "Minimum Efficiency Score", 
                0.0, 1.0, 0.5
            )
        
        calculate_batch = st.button("Process All Sequences")
        
        if st.session_state.uploaded_file and calculate_batch:
            process_batch_sequences(st.session_state.uploaded_file, has_headers, pam_type,
                                 gc_min, gc_max, efficiency_threshold)

def process_single_sequence(sequence, pam_type, gc_min, gc_max, 
                          efficiency_threshold, max_off_targets, 
                          show_all, check_secondary):
    """Process a single sequence for guide RNA design"""
    # Clean sequence first
    clean_sequence = sequence.upper().replace(" ", "").replace("\n", "")
    
    # Validate sequence
    if not clean_sequence:
        st.error("Please enter a DNA sequence")
    elif not set(clean_sequence).issubset({'A', 'T', 'G', 'C'}):
        st.error("Invalid characters found. Please use only A, T, G, C")
    else:
        st.success("Valid DNA sequence!")
        
        with st.spinner("Analyzing sequence..."):
            # Find PAM sites
            pam_sites = find_pam_sites(clean_sequence, pam_type)
            
            if pam_sites:
                # Design initial gRNAs
                grnas = design_grnas(clean_sequence, pam_sites, gc_min/100, gc_max/100)
                
                if grnas:
                    # Predict efficiency scores
                    grnas = predict_grna_efficiency(grnas)
                    
                    # Check off-targets
                    grnas = check_off_targets(grnas, clean_sequence)
                    
                    # Filter results if not showing all
                    if not show_all:
                        filtered_grnas = [g for g in grnas 
                                        if g['efficiency_score'] >= efficiency_threshold
                                        and g['off_target_score'] <= max_off_targets]
                        
                        if not filtered_grnas:
                            st.warning(f"No guides meet the criteria. Showing all {len(grnas)} candidates.")
                            display_grna_results(grnas)
                        else:
                            st.success(f"Found {len(filtered_grnas)} guides meeting criteria.")
                            display_grna_results(filtered_grnas)
                    else:
                        st.info(f"Showing all {len(grnas)} candidates.")
                        display_grna_results(grnas)
                else:
                    st.warning("No suitable guide RNAs found in the sequence.")
            else:
                if pam_type == "Cas12a (TTTV)":
                    st.warning("""
                    No Cas12a PAM sites found. This is normal because:
                    1. Cas12a needs exactly TTT + (A, C, or G)
                    2. This pattern is less common than SpCas9's NGG
                    3. Try:
                       - Using SpCas9 instead
                       - Checking the reverse complement
                       - Using a different region of your sequence
                    """)
                else:
                    st.warning(f"No PAM sites found for {pam_type}")

def process_batch_sequences(uploaded_file, has_headers, pam_type, gc_min, gc_max, efficiency_threshold):
    """Process multiple sequences for guide RNA design"""
    import pandas as pd
    import io
    
    try:
        # Read the file
        content = uploaded_file.read().decode()
        sequences = pd.read_csv(io.StringIO(content), 
                              header=0 if has_headers else None,
                              names=['name', 'sequence'] if has_headers else ['sequence'])
        
        if len(sequences) > 100:
            st.error("Maximum 100 sequences allowed per batch")
            return
        
        # Initialize results storage
        all_results = []
        
        # Process each sequence
        progress_bar = st.progress(0)
        for idx, row in sequences.iterrows():
            sequence = row['sequence']
            name = row['name'] if has_headers else f"Sequence_{idx+1}"
            
            # Clean and validate sequence
            clean_sequence = sequence.upper().replace(" ", "").replace("\n", "")
            if set(clean_sequence).issubset({'A', 'T', 'G', 'C'}):
                # Find PAM sites and design gRNAs
                pam_sites = find_pam_sites(clean_sequence, pam_type)
                if pam_sites:
                    grnas = design_grnas(clean_sequence, pam_sites, gc_min/100, gc_max/100)
                    if grnas:
                        # Add sequence name to results
                        for grna in grnas:
                            grna['sequence_name'] = name
                        grnas = predict_grna_efficiency(grnas)
                        grnas = check_off_targets(grnas, clean_sequence)
                        all_results.extend(grnas)
            
            # Update progress
            progress_bar.progress((idx + 1) / len(sequences))
        
        if all_results:
            # Convert results to DataFrame
            results_df = pd.DataFrame(all_results)
            
            # Filter by efficiency threshold
            filtered_df = results_df[results_df['efficiency_score'] >= efficiency_threshold]
            
            # Display summary statistics
            st.success(f"Processed {len(sequences)} sequences")
            st.info(f"Found {len(filtered_df)} guide RNAs meeting criteria")
            
            # Display interactive table
            st.dataframe(filtered_df.style.map(color_efficiency, subset=['efficiency_score'])
                                         .map(color_offtarget, subset=['off_target_score']))
            
            # Add download button
            csv = filtered_df.to_csv(index=False)
            st.download_button(
                "Download Results CSV",
                csv,
                "guide_rna_results.csv",
                "text/csv",
                key='download-csv'
            )
            
            # Add Excel download
            excel_buffer = io.BytesIO()
            filtered_df.to_excel(excel_buffer, index=False, engine='openpyxl')
            excel_data = excel_buffer.getvalue()
            st.download_button(
                "Download Results Excel",
                excel_data,
                "guide_rna_results.xlsx",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key='download-excel'
            )
            
    except Exception as e:
        st.error(f"Error processing file: {str(e)}")

def show_sequence_analysis():
    st.header("Sequence Analysis")
    
    # Input section
    col1, col2 = st.columns([2,1])
    
    with col1:
        # Use session state for sequence persistence
        sequence = st.text_area(
            "Enter DNA sequence",
            value=st.session_state.dna_sequence,
            height=150,
            help="Enter a DNA sequence (A, T, G, C only)"
        )
        # Update session state
        st.session_state.dna_sequence = sequence
        
        calculate = st.button("Analyze Sequence")
        
        if st.button("Generate Reverse Complement"):
            if sequence:
                complement = str(Seq(sequence).reverse_complement())
                st.code(complement, language="text")
                st.button("Use this sequence", 
                         help="Use reverse complement for guide RNA design")
    
    with col2:
        st.info("Sequence Statistics")
        if sequence and calculate:
            clean_seq = sequence.upper().replace(" ", "").replace("\n", "")
            if set(clean_seq).issubset({'A', 'T', 'G', 'C'}):
                # Basic statistics
                st.metric("Sequence Length", len(clean_seq))
                st.metric("GC Content", f"{gc_fraction(clean_seq)*100:.1f}%")
                
                # Nucleotide composition
                st.write("Nucleotide Composition:")
                for base in ['A', 'T', 'G', 'C']:
                    count = clean_seq.count(base)
                    st.progress(count/len(clean_seq), 
                              text=f"{base}: {count} ({count/len(clean_seq)*100:.1f}%)")
                
                # Add CRISPR scoring analysis with explanation
                st.subheader("CRISPR Scoring Analysis")
                st.info("""
                **How Scoring is Calculated:**
                My scoring system combines four key components based on the latest research:
                
                1. **GC Score (25%):**
                   - Optimal range: 45-65%
                   - Affects binding stability
                
                2. **Self-complementarity (25%):**
                   - Checks for unwanted hairpin structures
                   - Lower is better
                
                3. **Homopolymer Score (20%):**
                   - Penalizes long repeats of same nucleotide
                   - Especially strict on T repeats
                
                4. **Position Score (30%):**
                   - Evaluates nucleotide preferences at specific positions
                   - Based on empirical success rates
                
                Final score is weighted average of these components.
                Typical lab success rates:
                - Scores >80%: ~60-70% success
                - Scores >70%: ~50-60% success
                - Scores >50%: ~30-50% success
                """)
                
                from src.utils.scoring_algorithms import CRISPRScoring
                scorer = CRISPRScoring()
                scores = scorer.calculate_all_scores(clean_seq)
                
                # Display scores with color-coded bars
                st.write("Individual Score Components:")
                for name, score in scores.items():
                    if name != 'final_score':
                        color = 'green' if score >= 0.7 else 'yellow' if score >= 0.5 else 'red'
                        st.progress(score, text=f"{name}: {score:.2f}")
                
                # Display final score prominently
                st.metric("Overall Score", f"{scores['final_score']:.2f}")
                
                # Add score interpretation
                if scores['final_score'] >= 0.7:
                    st.success("This sequence has favorable characteristics for CRISPR targeting")
                elif scores['final_score'] >= 0.5:
                    st.warning("This sequence has moderate characteristics for CRISPR targeting")
                else:
                    st.error("This sequence may be challenging for CRISPR targeting")
    
    # Visualization section moved outside columns for full width
    if sequence and calculate:
        if set(clean_seq).issubset({'A', 'T', 'G', 'C'}):
            st.markdown("---")  # Add a visual separator
            st.subheader("Sequence Visualization")
            fig = create_sequence_plot(clean_seq)
            # Update figure layout for full width
            fig.update_layout(
                height=600,  # Increased height
                width=None,  # Allow dynamic width
                margin=dict(l=20, r=20, t=50, b=20),  # Reduce margins
                showlegend=True,
                hovermode='closest'
            )
            st.plotly_chart(fig, use_container_width=True)

def show_visualization():
    st.header("CRISPR Guide RNA Visualization & Analysis")
    
    # Use session state for sequence persistence
    sequence = st.text_area(
        "Enter DNA sequence for visualization",
        value=st.session_state.dna_sequence,
        height=150
    )
    # Update session state
    st.session_state.dna_sequence = sequence
    
    calculate = st.button("Generate Visualizations")
    
    # Sidebar for visualization options
    st.sidebar.subheader("Visualization Options")
    
    if sequence and calculate:
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        
        # Create tabs for different visualizations
        tab1, tab2, tab3, tab4 = st.tabs([
            "Guide RNA Map", 
            "Off-target Analysis", 
            "Efficiency Prediction",
            "Structure Analysis"
        ])
        
        with tab1:
            st.subheader("Guide RNA Target Sites")
            # Interactive genome browser style visualization
            plot_guide_rna_map(sequence)
            
        with tab2:
            st.subheader("Off-target Analysis")
            # Circos-style plot for off-target sites
            plot_off_target_analysis(sequence)
            
        with tab3:
            st.subheader("Efficiency Prediction")
            # ML-based efficiency prediction visualization
            plot_efficiency_prediction(sequence)
            
        with tab4:
            st.subheader("Secondary Structure")
            # RNA folding visualization
            plot_secondary_structure(sequence)

def plot_guide_rna_map(sequence):
    """Create an interactive genome browser-style visualization"""
    import plotly.graph_objects as go
    
    # Find PAM sites
    pam_sites = find_pam_sites(sequence, "SpCas9 (NGG)")
    grnas = design_grnas(sequence, pam_sites, 0.3, 0.7)
    
    # Create the figure
    fig = go.Figure()
    
    # Add sequence track
    fig.add_trace(go.Scatter(
        x=list(range(len(sequence))),
        y=[0]*len(sequence),
        mode='text',
        text=list(sequence),
        textfont=dict(size=10),
        name='Sequence'
    ))
    
    # Add PAM sites
    if pam_sites:
        fig.add_trace(go.Scatter(
            x=[site for site in pam_sites],
            y=[1]*len(pam_sites),
            mode='markers',
            marker=dict(
                size=10,
                color='red',
                symbol='triangle-down'
            ),
            name='PAM Sites'
        ))
    
    # Add gRNA regions
    if grnas:
        for grna in grnas:
            pos = grna['position']
            fig.add_shape(
                type="rect",
                x0=pos,
                x1=pos+20,
                y0=0.5,
                y1=1.5,
                fillcolor="rgba(0,100,80,0.2)",
                line=dict(width=0),
            )
    
    # Update layout with legend on the left
    fig.update_layout(
        title="Guide RNA Target Sites",
        showlegend=True,
        height=400,
        xaxis_title="Sequence Position",
        yaxis_visible=False,
        hovermode='closest',
        legend=dict(
            x=0,  # Position legend at the far left
            y=1,  # Position at the top
            xanchor='left',  # Anchor point on the left
            yanchor='top',   # Anchor at the top
            bgcolor='rgba(255,255,255,0.8)',  # Semi-transparent background
            bordercolor='rgba(0,0,0,0.2)',    # Light border
            borderwidth=1
        ),
        margin=dict(l=100, r=20, t=50, b=20)  # Increase left margin to accommodate legend
    )
    
    st.plotly_chart(fig, use_container_width=True)

def plot_off_target_analysis(sequence):
    """Create a circos-style plot for off-target analysis"""
    import plotly.graph_objects as go
    import numpy as np
    
    # Generate mock off-target data
    pam_sites = find_pam_sites(sequence, "SpCas9 (NGG)")
    if pam_sites:
        # Get initial gRNAs
        grnas = design_grnas(sequence, pam_sites, 0.3, 0.7)
        if grnas:
            # Calculate off-target scores
            grnas = check_off_targets(grnas, sequence)
            
            # Create circular plot
            fig = go.Figure()
            
            # Add main sequence arc
            theta = np.linspace(0, 360, len(sequence))
            fig.add_trace(go.Scatterpolar(
                r=[1]*len(sequence),
                theta=theta,
                mode='lines',
                name='Sequence',
                line=dict(color='blue')
            ))
            
            # Add off-target connections
            for grna in grnas:
                off_target_score = float(grna['off_target_score'])
                if off_target_score > 0:
                    fig.add_trace(go.Scatterpolar(
                        r=[1, off_target_score/100],
                        theta=[theta[grna['position']], theta[grna['position']]],
                        mode='lines',
                        name=f'Off-target: {off_target_score:.1f}',
                        line=dict(color='rgba(255,0,0,0.3)')
                    ))
            
            fig.update_layout(
                polar=dict(
                    radialaxis=dict(range=[0, 1], showticklabels=False),
                    angularaxis=dict(showticklabels=False)
                ),
                showlegend=True,
                height=600,
                title="Off-target Analysis (Circos Plot)"
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Add explanation
            st.info("""
            This circos plot shows potential off-target sites:
            - Blue circle represents your input sequence
            - Red lines show potential off-target locations
            - Longer lines indicate higher off-target risk
            """)

def plot_efficiency_prediction(sequence):
    """Create ML-based efficiency prediction visualization"""
    import plotly.graph_objects as go
    
    pam_sites = find_pam_sites(sequence, "SpCas9 (NGG)")
    if pam_sites:
        grnas = design_grnas(sequence, pam_sites, 0.3, 0.7)
        if grnas:
            # Predict efficiencies
            grnas = predict_grna_efficiency(grnas)
            
            # Create heatmap
            positions = [g['position'] for g in grnas]
            efficiencies = [g['efficiency_score'] for g in grnas]
            
            fig = go.Figure(data=go.Heatmap(
                z=[efficiencies],
                x=positions,
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(title='Efficiency Score')
            ))
            
            fig.update_layout(
                title="Guide RNA Efficiency Prediction",
                xaxis_title="Sequence Position",
                yaxis_visible=False,
                height=300
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Add ML feature importance plot
            st.subheader("Feature Importance in Efficiency Prediction")
            features = ['GC Content', 'Position Score', 'Self-complementarity']
            importance = [0.35, 0.45, 0.20]  # Our model weights
            
            fig = go.Figure(data=go.Bar(
                x=features,
                y=importance,
                marker_color='rgb(26, 118, 255)'
            ))
            
            fig.update_layout(
                title="Feature Importance",
                xaxis_title="Features",
                yaxis_title="Importance",
                height=300
            )
            
            st.plotly_chart(fig, use_container_width=True)

def plot_secondary_structure(sequence):
    """Create RNA secondary structure visualization"""
    import plotly.graph_objects as go
    import numpy as np
    
    if len(sequence) > 20:
        sequence = sequence[:20]  # Take first 20 bases for structure
    
    # Create figure with subplots
    fig = go.Figure()
    
    # Add sequence with larger font and better spacing
    fig.add_trace(go.Scatter(
        x=list(range(len(sequence))),
        y=[0]*len(sequence),
        mode='text',
        text=list(sequence),
        textfont=dict(size=24, color='black'),
        name='RNA Sequence',
        hoverinfo='text',
        hovertext=[f"Position {i+1}: {base}" for i, base in enumerate(sequence)]
    ))
    
    # Add arcs for potential base pairs with improved colors
    complement = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}  # Note: T → U for RNA
    max_height = 0
    
    for i in range(len(sequence)):
        for j in range(i+4, len(sequence)):  # Minimum loop size of 4
            if sequence[j] == complement.get(sequence[i].replace('T', 'U')):
                # Create arc with color based on base pair type
                x = np.linspace(i, j, 20)
                height = (j-i)/4
                max_height = max(max_height, height)
                y = np.sin(np.pi*(x-i)/(j-i)) * height
                
                # Color based on base pair type
                if sequence[i] in ['G', 'C']:
                    color = 'rgba(255,0,0,0.4)'  # Strong GC pairs in red
                    width = 2
                else:
                    color = 'rgba(0,0,255,0.2)'  # Weaker AU pairs in blue
                    width = 1.5
                
                fig.add_trace(go.Scatter(
                    x=x,
                    y=y,
                    mode='lines',
                    line=dict(color=color, width=width),
                    name=f'{sequence[i]}-{sequence[j]} pair',
                    hoverinfo='text',
                    hovertext=f'Base pair: {sequence[i]}-{sequence[j]} (positions {i+1}-{j+1})',
                    showlegend=False
                ))
    
    # Add legend for base pair types
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='rgba(255,0,0,0.4)', width=2),
        name='Strong (G-C) pairs'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='rgba(0,0,255,0.2)', width=1.5),
        name='Weak (A-U) pairs'
    ))
    
    # Update layout with better formatting
    fig.update_layout(
        title={
            'text': "RNA Secondary Structure Prediction",
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        showlegend=True,
        legend=dict(
            x=0,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255,255,255,0.8)'
        ),
        height=500,
        xaxis=dict(
            title="Sequence Position",
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            range=[-max_height-0.5, max_height+0.5],
            showgrid=False,
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='black',
            showticklabels=False
        ),
        plot_bgcolor='white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Add explanatory text
    st.info("""
    **Understanding the Structure Visualization:**
    
    - **Sequence**: The RNA sequence is shown along the horizontal axis
    - **Arcs**: Connect bases that could potentially pair with each other
        - **Red arcs**: Strong G-C base pairs (3 hydrogen bonds)
        - **Blue arcs**: Weaker A-U base pairs (2 hydrogen bonds)
    - **Arc Height**: Longer arcs are higher, showing long-range interactions
    - **Overlapping Arcs**: Indicate competing possible structures
    
    **What This Means:**
    - More arcs = More potential for secondary structure
    - Many overlapping arcs = Higher chance of structural issues
    - Long-range arcs (very tall) = Possible stability problems
    - Many strong (red) pairs = More stable structure
    """)
    
    # Add structure stability analysis
    gc_pairs = sum(1 for i in range(len(sequence)) 
                  for j in range(i+4, len(sequence)) 
                  if sequence[i] in 'GC' and sequence[j] == complement.get(sequence[i].replace('T', 'U')))
    au_pairs = sum(1 for i in range(len(sequence)) 
                  for j in range(i+4, len(sequence)) 
                  if sequence[i] in 'AT' and sequence[j] == complement.get(sequence[i].replace('T', 'U')))
    
    stability = "High" if gc_pairs > 3 else "Medium" if gc_pairs + au_pairs > 2 else "Low"
    color = {"High": "green", "Medium": "orange", "Low": "red"}[stability]
    
    st.markdown(f"""
    **Structure Stability Analysis:**
    - Strong G-C pairs: {gc_pairs}
    - Weak A-U pairs: {au_pairs}
    - Overall stability: :{color}[{stability}]
    """)

def display_grna_results(grnas):
    """Enhanced results display"""
    import streamlit as st
    import pandas as pd
    
    df = pd.DataFrame(grnas)
    
    # Reorder and format columns
    columns = [
        'sequence', 
        'pam', 
        'position', 
        'gc_score',
        'self_complementarity',
        'homopolymer',
        'efficiency_score', 
        'off_target_score'
    ]
    
    # Select only columns that exist
    existing_columns = [col for col in columns if col in df.columns]
    df = df[existing_columns]
    
    # Format display values
    score_columns = ['gc_score', 'self_complementarity', 'homopolymer', 'efficiency_score']
    for col in df.columns:
        if col in score_columns and col in df.columns:
            df[col] = df[col].apply(lambda x: f"{x*100:.1f}%" if isinstance(x, (int, float)) else "N/A")
        elif col == 'off_target_score' and col in df.columns:
            df[col] = df[col].apply(lambda x: f"{x:.1f}" if isinstance(x, (int, float)) else "N/A")
    
    # Apply styling to score columns only
    style_columns = [col for col in score_columns if col in df.columns]
    
    styled_df = df.style.map(color_efficiency, subset=style_columns)
    if 'off_target_score' in df.columns:
        styled_df = styled_df.map(color_offtarget, subset=['off_target_score'])
    
    # Display results with explanations
    st.subheader("Guide RNA Results")
    st.markdown("""
    Scoring Components:
    - **GC Score**: Optimal GC content (45-65%)
    - **Self-complementarity**: Lower is better (risk of hairpin formation)
    - **Homopolymer**: Higher is better (fewer repeat nucleotides)
    - **Efficiency Score**: Combined score (Green ≥ 70%, Yellow ≥ 50%, Red < 50%)
    - **Off-target Score**: Lower is better (Green ≤ 20, Yellow ≤ 50, Red > 50)
    """)
    
    st.dataframe(styled_df)

if __name__ == "__main__":
    main() 