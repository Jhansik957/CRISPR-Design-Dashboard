import re
from Bio.SeqUtils import gc_fraction
import numpy as np
from sklearn.preprocessing import StandardScaler
import streamlit as st
from .scoring_algorithms import CRISPRScoring

def find_pam_sites(sequence, pam_type):
    """Improved PAM site finding with better sequence handling"""
    # Clean the sequence first
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    
    # Updated PAM patterns
    pam_patterns = {
        "SpCas9 (NGG)": "(?=.GG)",
        "SaCas9 (NNGRRT)": "(?=..[GA][AG][AG]T)",
        "Cas12a (TTTV)": "(?=TTT[ACG])"  # More precise TTTV pattern
    }
    
    pattern = pam_patterns[pam_type]
    sites = [m.start() for m in re.finditer(pattern, sequence)]
    
    # Add debugging information
    if not sites:
        if pam_type == "Cas12a (TTTV)":
            st.info("""
            No Cas12a PAM sites (TTTV) found. This is common because:
            - Needs exactly TTT followed by A, C, or G
            - More restrictive than SpCas9 or SaCas9
            - Try using the reverse complement of your sequence
            """)
    
    return sites

def design_grnas(sequence, pam_sites, gc_min, gc_max):
    grnas = []
    for site in pam_sites:
        # For SpCas9, gRNA is 20nt upstream of PAM
        if site >= 20:
            grna_seq = sequence[site-20:site]
            gc_content = gc_fraction(grna_seq)
            
            if gc_min <= gc_content <= gc_max:
                grnas.append({
                    'sequence': grna_seq,
                    'pam': sequence[site:site+3],
                    'position': site-20,
                    'gc_content': gc_content
                })
    return grnas

def analyze_secondary_structure(grnas):
    """Simple secondary structure analysis without Vienna RNA"""
    for grna in grnas:
        sequence = grna['sequence']
        
        # Simple estimation of structure stability
        stability_score = simple_structure_estimation(sequence)
        
        grna['secondary_structure'] = "Estimated"
        grna['free_energy'] = stability_score
    
    return grnas

def simple_structure_estimation(sequence):
    """Simplified stability estimation based on base pairing potential"""
    # Count potential G-C pairs (stronger bonds)
    gc_pairs = min(sequence.count('G'), sequence.count('C'))
    
    # Count potential A-T pairs (weaker bonds)
    at_pairs = min(sequence.count('A'), sequence.count('T'))
    
    # Simple scoring: GC pairs contribute more to stability
    stability = -(gc_pairs * 3 + at_pairs * 2)
    return round(stability, 2)

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
        'gc_content',
        'efficiency_score', 
        'off_target_score'
    ]
    
    df = df[columns]
    
    # Format display values
    df['gc_content'] = df['gc_content'].apply(lambda x: f"{x*100:.1f}%")
    df['efficiency_score'] = df['efficiency_score'].apply(lambda x: f"{x*100:.1f}%")
    df['off_target_score'] = df['off_target_score'].apply(lambda x: f"{x:.1f}")
    
    # Add color coding
    def color_efficiency(val):
        try:
            num_val = float(val.strip('%'))/100
            if num_val >= 0.7:
                return 'background-color: #90EE90'  # Light green
            elif num_val >= 0.5:
                return 'background-color: #FFFFE0'  # Light yellow
            return 'background-color: #FFB6C1'      # Light red
        except:
            return ''

    def color_offtarget(val):
        try:
            num_val = float(val)
            if num_val <= 20:
                return 'background-color: #90EE90'  # Light green
            elif num_val <= 50:
                return 'background-color: #FFFFE0'  # Light yellow
            return 'background-color: #FFB6C1'      # Light red
        except:
            return ''
    
    # Apply styling
    styled_df = df.style.map(color_efficiency, subset=['efficiency_score'])\
                       .map(color_offtarget, subset=['off_target_score'])
    
    # Display results with explanations
    st.subheader("Guide RNA Results")
    st.markdown("""
    - **Efficiency Score**: Higher is better (Green ≥ 70%, Yellow ≥ 50%, Red < 50%)
    - **Off-target Score**: Lower is better (Green ≤ 20, Yellow ≤ 50, Red > 50)
    """)
    
    st.dataframe(styled_df)

def predict_grna_efficiency(grnas):
    """Predict gRNA efficiency using our custom scoring algorithm"""
    from .scoring_algorithms import CRISPRScoring
    scorer = CRISPRScoring()
    
    for grna in grnas:
        sequence = grna['sequence']
        # Calculate all scores using our new system
        scores = scorer.calculate_all_scores(sequence)
        # Use the final_score as our efficiency score
        grna['efficiency_score'] = scores['final_score']
        # Store individual component scores for reference
        grna['gc_score'] = scores['gc_score']
        grna['self_complementarity'] = scores['self_complementarity']
        grna['homopolymer'] = scores['homopolymer']
        grna['position'] = scores['position']
    
    return grnas

def analyze_structure_details(sequence):
    """Detailed structure analysis"""
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    
    return {
        'gc_pairs': gc_count,
        'at_pairs': at_count,
        'stability': 'High' if gc_count > len(sequence)/2 else 'Medium' if gc_count > len(sequence)/3 else 'Low'
    }

def calculate_position_score(sequence):
    """Calculate position-specific nucleotide preferences"""
    # Position-specific weights based on literature
    weights = {
        'G': [0.8, 0.6, 0.7, 0.5, 0.6, 0.7, 0.8, 0.9, 0.8, 0.7,
              0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5],
        'A': [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.7, 0.6, 0.5,
              0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        'T': [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4,
              0.5, 0.6, 0.7, 0.8, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4],
        'C': [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.8, 0.7, 0.6,
              0.5, 0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    }
    
    score = 0
    for i, base in enumerate(sequence):
        score += weights[base][i]
    
    return score / len(sequence)

def check_self_complementarity(sequence):
    """Check for self-complementarity issues"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = ''.join(complement[base] for base in reversed(sequence))
    
    max_complementary = 0
    for i in range(len(sequence)-3):
        for j in range(i+3, len(sequence)):
            if sequence[i:j] in rev_comp:
                max_complementary = max(max_complementary, j-i)
    
    return max_complementary / len(sequence)

def check_off_targets(grnas, target_sequence):
    """Improved off-target prediction"""
    for grna in grnas:
        sequence = grna['sequence']
        
        # Initialize off-target score
        off_target_score = 0
        
        # Scan through sequence with sliding window
        for i in range(len(target_sequence)-len(sequence)):
            window = target_sequence[i:i+len(sequence)]
            mismatches = sum(a != b for a, b in zip(sequence, window))
            
            # Only consider windows with 4 or fewer mismatches
            if mismatches <= 4:
                # Weight score by number of mismatches
                # Fewer mismatches = higher off-target risk
                off_target_score += 1.0 / (2 ** mismatches)
        
        # Normalize score to 0-100 range
        normalized_score = min(100, off_target_score * 20)
        grna['off_target_score'] = round(normalized_score, 1)
    
    return grnas 

def create_sequence_plot(sequence):
    """Create an interactive sequence visualization"""
    import plotly.graph_objects as go
    
    # Create base colors
    colors = {
        'A': '#FF9999', 'T': '#99FF99',
        'G': '#9999FF', 'C': '#FFFF99'
    }
    
    # Create hover text
    hover_text = [f"Position: {i+1}<br>Base: {base}" 
                 for i, base in enumerate(sequence)]
    
    # Create the figure
    fig = go.Figure()
    
    # Add bases as colored points
    for base in 'ATGC':
        positions = [i for i, b in enumerate(sequence) if b == base]
        if positions:
            fig.add_trace(go.Scatter(
                x=positions,
                y=[1]*len(positions),
                mode='markers',
                name=base,
                marker=dict(
                    size=20,
                    color=colors[base],
                    line=dict(width=1, color='black')
                ),
                text=[hover_text[p] for p in positions],
                hoverinfo='text'
            ))
    
    # Update layout
    fig.update_layout(
        title="DNA Sequence Visualization",
        xaxis_title="Position",
        yaxis_visible=False,
        height=200,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    return fig 