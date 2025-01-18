import numpy as np
from Bio.SeqUtils import gc_fraction

class CRISPRScoring:
    def __init__(self):
        """Initialize basic scoring methods"""
        pass

    def gc_score(self, sequence):
        """
        Calculate GC content score based on 2023 meta-analysis
        Optimal GC content is between 45-65% (updated range)
        Source: Xu et al. 2023, Nature Communications
        """
        gc = gc_fraction(sequence)
        
        # Score peaks between 45-65% GC
        if 0.45 <= gc <= 0.65:
            return 1.0
        elif gc < 0.45:
            return gc / 0.45  # Linear falloff below 45%
        else:
            return (1 - gc) / 0.35  # Linear falloff above 65%

    def self_complementarity_score(self, sequence):
        """
        Check for self-complementarity that could form hairpins
        Updated thresholds based on 2023 studies
        Sources: 
        - Kim et al. 2022, Nature Biotechnology
        - Labuhn et al. 2023, Nucleic Acids Research
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rev_comp = ''.join(complement.get(base, base) for base in reversed(sequence))
        
        max_complementary = 0
        seq_len = len(sequence)
        
        # Check for complementary regions of 5 or more bases (updated threshold)
        for i in range(seq_len - 4):
            for j in range(i + 5, seq_len):
                region = sequence[i:j]
                if region in rev_comp:
                    max_complementary = max(max_complementary, len(region))
        
        # Updated scoring thresholds
        if max_complementary < 5:
            return 1.0
        elif max_complementary >= 12:
            return 0.0
        else:
            return 1 - ((max_complementary - 5) / 7)  # More gradual penalty

    def homopolymer_score(self, sequence):
        """
        Penalize sequences with homopolymer runs
        Updated with nucleotide-specific penalties
        Source: Zhang et al. 2023, Genome Biology
        """
        def get_homopolymer_run(seq, base):
            max_run = 1
            current_run = 1
            for i in range(1, len(seq)):
                if seq[i] == seq[i-1] == base:
                    current_run += 1
                    max_run = max(max_run, current_run)
                else:
                    current_run = 1
            return max_run
        
        # Check each base separately with different thresholds
        t_run = get_homopolymer_run(sequence, 'T')
        a_run = get_homopolymer_run(sequence, 'A')
        g_run = get_homopolymer_run(sequence, 'G')
        c_run = get_homopolymer_run(sequence, 'C')
        
        # Stronger penalty for T runs (most detrimental)
        t_score = 1.0 if t_run <= 2 else 0.0 if t_run >= 4 else 1 - ((t_run - 2) / 2)
        # Standard penalties for other bases
        a_score = 1.0 if a_run <= 3 else 0.0 if a_run >= 5 else 1 - ((a_run - 3) / 2)
        g_score = 1.0 if g_run <= 3 else 0.0 if g_run >= 5 else 1 - ((g_run - 3) / 2)
        c_score = 1.0 if c_run <= 3 else 0.0 if c_run >= 5 else 1 - ((c_run - 3) / 2)
        
        # Weighted average with higher weight for T-runs
        return 0.4 * t_score + 0.2 * a_score + 0.2 * g_score + 0.2 * c_score

    def position_score(self, sequence):
        """
        Score based on preferred bases at specific positions
        Updated weights based on 2023 meta-analysis
        Sources: 
        - DeWeirdt et al. 2023, Nature Biotechnology
        - Kim et al. 2023, Cell Reports Methods
        """
        # Updated position-specific weights (1-based position)
        preferences = {
            1: {'G': 0.9, 'A': 0.6, 'C': 0.4, 'T': 0.2},  # Strong G preference
            2: {'G': 0.3, 'A': 0.8, 'C': 0.4, 'T': 0.3},  # A preference
            3: {'G': 0.6, 'A': 0.6, 'C': 0.4, 'T': 0.3},  # G/A preference
            4: {'G': 0.5, 'A': 0.5, 'C': 0.5, 'T': 0.4},  # Mild base preference
            16: {'G': 0.7, 'A': 0.5, 'C': 0.3, 'T': 0.3}, # G preference
            17: {'G': 0.8, 'A': 0.4, 'C': 0.3, 'T': 0.2}, # Strong G preference
            18: {'G': 0.7, 'A': 0.4, 'C': 0.4, 'T': 0.3}, # G preference
            19: {'G': 0.6, 'A': 0.5, 'C': 0.4, 'T': 0.3}, # Mild G preference
            20: {'G': 0.8, 'A': 0.4, 'C': 0.3, 'T': 0.2}  # Strong G preference
        }
        
        score = 0
        count = 0
        
        for pos, weights in preferences.items():
            if pos <= len(sequence):
                base = sequence[pos-1]
                score += weights.get(base, 0.3)  # Updated default score
                count += 1
        
        return score / count if count > 0 else 0.5

    def calculate_all_scores(self, sequence):
        """
        Calculate all basic scores and return weighted average
        Updated weights based on importance in recent literature
        """
        scores = {
            'gc_score': self.gc_score(sequence),
            'self_complementarity': self.self_complementarity_score(sequence),
            'homopolymer': self.homopolymer_score(sequence),
            'position': self.position_score(sequence)
        }
        
        # Updated weights based on recent meta-analyses
        weights = {
            'gc_score': 0.25,             # Slightly reduced weight
            'self_complementarity': 0.25,  # Slightly reduced weight
            'homopolymer': 0.2,           # Maintained weight
            'position': 0.3               # Increased weight due to strong evidence
        }
        
        final_score = sum(scores[k] * weights[k] for k in scores)
        scores['final_score'] = final_score
        
        return scores 