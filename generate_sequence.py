import random
import streamlit as st

def generate_dna_sequence(length, organism_gc):
    """Generate a random DNA sequence with biologically relevant GC content"""
    gc_bias = organism_gc / 100  # Convert percentage to decimal
    
    # Weights based on biological GC content
    weights = {
        'A': (1 - gc_bias) / 2,
        'T': (1 - gc_bias) / 2,
        'G': gc_bias / 2,
        'C': gc_bias / 2
    }
    
    nucleotides = list(weights.keys())
    probabilities = list(weights.values())
    
    sequence = ''.join(random.choices(nucleotides, weights=probabilities, k=length))
    return sequence

def main():
    st.title("Biological DNA Sequence Generator")
    
    # Simple input for sequence length
    length = st.text_input(
        "Enter number of nucleotides",
        value="100"
    )
    
    # Organism-based GC content selection
    organism_type = st.selectbox(
        "Select Organism Type",
        options=[
            "Human (42% GC)",
            "E. coli (51% GC)",
            "S. cerevisiae (38% GC)",
            "Mycobacterium (65% GC)",
            "Custom GC%"
        ]
    )
    
    # GC content mapping
    gc_content_map = {
        "Human (42% GC)": 42,
        "E. coli (51% GC)": 51,
        "S. cerevisiae (38% GC)": 38,
        "Mycobacterium (65% GC)": 65,
        "Custom GC%": None
    }
    
    # Show custom slider only if Custom GC% is selected
    gc_percentage = gc_content_map[organism_type]
    if organism_type == "Custom GC%":
        gc_percentage = st.slider(
            "Custom GC Content (%)",
            min_value=25,
            max_value=75,
            value=50,
            help="Natural GC content ranges from 25-75% in different organisms"
        )
    
    if st.button("Generate Sequence"):
        try:
            length = int(length)
            if length <= 0:
                st.error("Please enter a positive number")
                return
                
            sequence = generate_dna_sequence(length, gc_percentage)
            
            # Display sequence
            st.subheader("Generated DNA Sequence:")
            st.text_area("Sequence", sequence, height=150)
            
            # Detailed statistics in columns
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Length", len(sequence))
            with col2:
                gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
                st.metric("GC Content", f"{gc_content*100:.1f}%")
            with col3:
                st.metric("AT Content", f"{(1-gc_content)*100:.1f}%")
            with col4:
                st.metric("Nucleotides", 
                         f"A:{sequence.count('A')} T:{sequence.count('T')} "
                         f"G:{sequence.count('G')} C:{sequence.count('C')}")
            
            # Add biological context
            st.info(f"""
            This sequence has characteristics similar to {organism_type.split('(')[0]} DNA:
            - Target GC content: {gc_percentage}%
            - Actual GC content: {gc_content*100:.1f}%
            - Natural range: Most organisms have 35-55% GC content
            """)
            
            # Improved copy functionality
            st.text("Click below to copy the sequence:")
            st.code(sequence, language="text")  # Display in a copyable code block
            
            # Add download button
            st.download_button(
                label="Download Sequence",
                data=sequence,
                file_name="dna_sequence.txt",
                mime="text/plain"
            )

        except ValueError:
            st.error("Please enter a valid number")

if __name__ == "__main__":
    main() 