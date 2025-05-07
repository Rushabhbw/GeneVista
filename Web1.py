import streamlit as st
from Bio import Entrez, SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Data import CodonTable
from io import StringIO
import pandas as pd
import base64
from streamlit_extras.stylable_container import stylable_container
from PIL import Image
import requests
import time
import re

# ====== GENETIC CODE MAPPING ======
GENETIC_CODE_NAMES = {
    1: "Standard",
    2: "Vertebrate Mitochondrial", 
    3: "Yeast Mitochondrial",
    4: "Mold/Protozoan Mitochondrial",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate Nuclear",
    9: "Echinoderm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Flatworm Mitochondrial",
    16: "Chlorophycean Mitochondrial",
    21: "Trematode Mitochondrial",
    22: "Scenedesmus Mitochondrial",
    23: "Thraustochytrium Mitochondrial"
}

# Set your email for NCBI Entrez
Entrez.email = "wakaderushabh@gmail.com"

# Streamlit page configuration
st.set_page_config(page_title="🔬 GeneVista: Bioinformatics Sequence Hub", layout="wide")

# ====== ENHANCED VISUAL DESIGN ======
st.markdown(f"""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Poppins:wght@300;400;600;700&display=swap');
    
    /* Main app styling */
    html, body, [class*="css"] {{
        font-family: 'Poppins', sans-serif;
    }}
    
    /* Animated gradient background */
    .stApp {{
        background: linear-gradient(-45deg, #f5f7fa, #e4e8f0, #f0f4f9, #e8ecf2);
        background-size: 400% 400%;
        animation: gradient 15s ease infinite;
    }}
    
    @keyframes gradient {{
        0% {{ background-position: 0% 50%; }}
        50% {{ background-position: 100% 50%; }}
        100% {{ background-position: 0% 50%; }}
    }}
    
    /* Content cards */
    .main-content {{
        background: rgba(255, 255, 255, 0.93);
        border-radius: 15px;
        padding: 2rem;
        margin: 1rem 0;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        border: 1px solid rgba(255, 255, 255, 0.2);
        animation: fadeIn 1s ease-out;
    }}
    
    /* Navigation buttons */
    .stButton>button {{
        border-radius: 12px !important;
        padding: 10px 24px !important;
        transition: all 0.3s !important;
        background: rgba(255,255,255,0.9) !important;
        border: 1px solid #e0e0e0 !important;
    }}
    
    .stButton>button:hover {{
        transform: translateY(-2px) !important;
        box-shadow: 0 6px 12px rgba(0,0,0,0.1) !important;
    }}
    
    /* Input fields */
    .stTextInput>div>div>input, 
    .stTextArea>textarea {{
        background: rgba(255,255,255,0.95) !important;
        border-radius: 12px !important;
        padding: 12px !important;
    }}
    
    /* Custom footer */
    .custom-footer {{
        margin-top: 50px;
        padding: 15px;
        color: white;
        background: rgba(44, 62, 80, 0.9);
        border-radius: 12px;
        backdrop-filter: blur(5px);
    }}
    
    /* Animation for content */
    @keyframes fadeIn {{
        0% {{ opacity: 0; transform: translateY(10px); }}
        100% {{ opacity: 1; transform: translateY(0); }}
    }}
    
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 8px;
    }}
    
    .stTabs [data-baseweb="tab"] {{
        background: rgba(255,255,255,0.8);
        border-radius: 8px 8px 0 0;
        padding: 10px 20px;
        transition: all 0.3s;
    }}
    
    .stTabs [aria-selected="true"] {{
        background: white;
        color: #4CAF50;
        font-weight: 600;
    }}
</style>
""", unsafe_allow_html=True)

# JavaScript for smooth transitions
st.markdown("""
<script>
document.addEventListener('DOMContentLoaded', () => {
    // Smooth transitions between pages
    const observer = new MutationObserver(() => {
        document.querySelectorAll('.main-content').forEach(el => {
            el.style.animation = 'fadeIn 0.8s ease-out';
        });
    });
    observer.observe(document.body, { childList: true, subtree: true });
});
</script>
""", unsafe_allow_html=True)

# --- Top Navigation State Handling ---
if "top_nav" not in st.session_state:
    st.session_state.top_nav = "🏠 Home"

# Render Top Navigation as Buttons
with st.container():
    col1, col2, col3, col4, col5, col6, col7 = st.columns([2, 1, 1, 1, 1, 1, 1])
    col1.markdown("<h1 style='margin-bottom:0;color:#2c3e50;'>GeneVista</h1>", unsafe_allow_html=True)
    if col2.button("🏠 Home"):
        st.session_state.top_nav = "🏠 Home"
    if col3.button("📖 User Guide"):
        st.session_state.top_nav = "📖 User Guide"
    if col4.button("🧬 Fetch Sequence"):
        st.session_state.top_nav = "🧬 Fetch Sequence (NCBI)"
    if col5.button("🔗 Alignment"):
        st.session_state.top_nav = "🔗 Sequence Alignment"
    if col6.button("🧬 Translation"):
        st.session_state.top_nav = "🧬 Sequence Translation"
    if col7.button("📌 About"):
        st.session_state.top_nav = "📌 About"

# Use top_nav state to control page
page = st.session_state.top_nav

# Optional: Hide Sidebar if not needed
st.markdown("""
    <style>
    [data-testid="stSidebar"] {
        display: none;
    }
    </style>
""", unsafe_allow_html=True)

# Function to create a download link for GenBank file or alignment result
def download_link(content, filename, label):
    b64 = base64.b64encode(content.encode()).decode()
    href = f'<a href="data:file/txt;base64,{b64}" download="{filename}">{label}</a>'
    return href

# ------------------------------ HOME ------------------------------
if page == "🏠 Home":
    with stylable_container(key="home_container", css_styles="""
        { background: rgba(255,255,255,0.93); border-radius: 15px; padding: 2rem; }
    """):
        st.title("🧬 GeneVista: Bioinformatics Sequence Hub")

        st.markdown("""
        Welcome to **GeneVista**, your integrated bioinformatics web server for essential sequence analysis tasks.  
        This platform is designed for **students, researchers, and educators** working in the fields of genomics, molecular biology, and computational biology.

        ---

        ### 🚀 Key Features:

        - **🔍 Sequence Retrieval from NCBI:**  
        Easily fetch nucleotide and protein sequences using accession numbers. View detailed metadata, taxonomy, GenBank records, and sequence statistics.

        - **🧩 Pairwise Sequence Alignment:**  
        Perform both **Global (Needleman-Wunsch)** and **Local (Smith-Waterman)** alignments for nucleotide or protein sequences. Customize scoring matrices, gap penalties, and visualize alignment metrics.

        - **🧬 Sequence Translation:**  
        Convert DNA/RNA to protein sequences using standard or specialized genetic codes (including mitochondrial and bacterial codes). View all 6 reading frames.

        - **📑 Downloadable Reports:**  
        Export alignment results, GenBank files, and protein translations for offline analysis and documentation.

        - **📖 Easy-to-Follow User Guide:**  
        A step-by-step manual to help you navigate and use the platform efficiently.

        ---

        ### 🎯 Why Use GeneVista?

        - **No installation needed** — fully browser-based.
        - Beginner-friendly, with expandable sections and interactive controls.
        - Direct access to public biological databases.
        - Designed for learning, research, and practical applications in bioinformatics.

        ---

        💡 *Explore the tools via the top navigation and start analyzing your sequences now!*
        """)

# ------------------------------ USER GUIDE ------------------------------
elif page == "📖 User Guide":
    with stylable_container(key="guide_container", css_styles="""
        { background: rgba(255,255,255,0.93); border-radius: 15px; padding: 2rem; }
    """):
        st.title("📖 User Guide")
        st.markdown("""
        ### 📌 Features Overview:

        1. **🧬 Fetch Sequence (NCBI)**
           - Input either nucleotide or protein accession number.
           - Choose database: `nucleotide` or `protein`.
           - View metadata, description, organism, sequence, and GenBank data.

        2. **🔗 Sequence Alignment**
           - Input sequences manually or upload FASTA files
           - Choose between Global (Needleman-Wunsch) and Local (Smith-Waterman) alignment
           - Choose sequence type: Nucleotide or Protein
           - Customize scoring and gap penalties
           - View aligned results including score, identity, and match statistics

        3. **🧬 Sequence Translation**
           - Input DNA/RNA sequence
           - Choose between standard translation or all 6 reading frames
           - Select from 23 NCBI genetic code tables
           - View protein sequence with stop codons marked (*)

        ### 🔗 Additional Resources:
        - [Biopython Documentation](https://biopython.org/wiki/Documentation)
        - [NCBI Handbook](https://www.ncbi.nlm.nih.gov/books/NBK143764/)
        - [Genetic Code Tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
        """)

# ------------------------------ FETCH SEQUENCE ------------------------------
elif page == "🧬 Fetch Sequence (NCBI)":
    with stylable_container(key="fetch_container", css_styles="""
        { background: rgba(255,255,255,0.93); border-radius: 15px; padding: 2rem; }
    """):
        st.title("🧬 Fetch Sequence from NCBI")
        st.markdown("Input a valid accession number (e.g., NM_001200001 for nucleotide or NP_001191970 for protein).")

        accession = st.text_input("🔑 Enter Accession Number")
        db_choice = st.selectbox("📚 Choose Database", ["nucleotide", "protein"])

        if st.button("🔍 Fetch Sequence") and accession:
            try:
                with Entrez.efetch(db=db_choice, id=accession, rettype="gb", retmode="text") as handle:
                    record = SeqIO.read(handle, "genbank")

                    st.markdown("## 🧬 Basic Metadata")
                    st.write(f"**Definition:** {record.description}")
                    st.write(f"**Organism:** {record.annotations.get('organism', 'N/A')}")
                    st.write(f"**Accession:** {accession}")
                    st.write(f"**Date:** {record.annotations.get('date', 'N/A')}")
                    st.write(f"**NCBI Link:** [View in NCBI](https://www.ncbi.nlm.nih.gov/nuccore/{accession})")

                    taxonomy = ", ".join(record.annotations.get("taxonomy", []))
                    st.write(f"**Taxonomy:** {taxonomy if taxonomy else 'N/A'}")
                    st.write(f"**Molecule Type:** {record.annotations.get('molecule_type', 'N/A')}")
                    st.write(f"**Topology:** {record.annotations.get('topology', 'N/A')}")

                    st.markdown("## 📊 Sequence Stats")
                    sequence = str(record.seq)
                    length = len(sequence)
                    gc_content = 100 * (sequence.count("G") + sequence.count("C")) / length if length > 0 else 0
                    st.write(f"**Length:** {length} bp")
                    st.write(f"**GC Content:** {gc_content:.2f}%")

                    with st.expander("🧬 Show Sequence (First 100 bp)"):
                        st.code(sequence[:100] + "...", language="text")

                    with st.expander("📖 Show Full Sequence"):
                        st.text_area("Full Sequence", sequence, height=300)

                    st.markdown("## 🔍 Features Table")
                    if record.features:
                        feature_data = []
                        for feature in record.features:
                            feature_data.append({
                                "Type": feature.type,
                                "Location": str(feature.location),
                                "Qualifiers": str(feature.qualifiers)
                            })
                        df = pd.DataFrame(feature_data)
                        st.dataframe(df)
                    else:
                        st.write("No features available.")

                    st.markdown("## 📄 Full GenBank Record")
                    output = StringIO()
                    SeqIO.write(record, output, "genbank")
                    genbank_text = output.getvalue()
                    st.text_area("🧬 GenBank File", genbank_text, height=300)
                    st.markdown(download_link(genbank_text, f"{accession}.gb", "📂 Download GenBank File"), unsafe_allow_html=True)

                    st.markdown("## 🧠 Gene/Protein Insights")
                    st.markdown("### 🌍 NCBI Summary")
                    st.markdown(f"🔗 [NCBI Search for {accession}](https://www.ncbi.nlm.nih.gov/search/all/?term={accession})")

                    st.markdown("### 🧪 UniProt Search")
                    st.markdown(f"🔗 [UniProt Entry for {accession}](https://www.uniprot.org/uniprotkb?query={accession})")

                    st.markdown("### 📚 PubMed Literature")
                    st.markdown(f"🔗 [PubMed Articles on {accession}](https://pubmed.ncbi.nlm.nih.gov/?term={accession})")

            except Exception as e:
                st.error(f"Error fetching data: {e}")

# ------------------------------ ALIGNMENT ------------------------------
elif page == "🔗 Sequence Alignment":
    with stylable_container(key="alignment_container", css_styles="""
        { background: rgba(255,255,255,0.93); border-radius: 15px; padding: 2rem; }
    """):
        st.title("🔗 Sequence Alignment")
        alignment_type = st.radio("Choose Alignment Type", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"])
        seq_type = st.radio("Select Sequence Type", ["Nucleotide", "Protein"])

        # New tab interface for input method
        tab1, tab2 = st.tabs(["✏️ Manual Entry", "📁 Upload FASTA Files"])
        
        seq1 = ""
        seq2 = ""
        
        with tab1:
            if seq_type == "Nucleotide":
                if st.button("🔁 Use Example Nucleotide Sequences"):
                    st.session_state.seq1 = "ATGCGTACGTTAG"
                    st.session_state.seq2 = "ATGCGTTCGTAG"
            else:
                if st.button("🔁 Use Example Protein Sequences"):
                    st.session_state.seq1 = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGQKL"
                    st.session_state.seq2 = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDLVTVEAAETFSLNNLGQKL"

            seq1 = st.text_area("🧬 Enter First Sequence", value=st.session_state.get("seq1", ""))
            seq2 = st.text_area("🧬 Enter Second Sequence", value=st.session_state.get("seq2", ""))

        with tab2:
            st.markdown("""
            **📌 Upload Options:**
            - Upload **two sequences in one FASTA file** (for pairwise alignment)
            - **OR** upload two separate FASTA files (one sequence per file)
            - Maximum file size: 5MB
            - Supported formats: `.fasta`, `.fa`, `.txt`
            """)
            
            col1, col2 = st.columns(2)
            with col1:
                combined_file = st.file_uploader("Combined FASTA (2 sequences)", 
                                              type=["fasta", "fa", "txt"],
                                              accept_multiple_files=False,
                                              help="Single file containing both sequences")
            with col2:
                separate_file1 = st.file_uploader("First Sequence FASTA", 
                                               type=["fasta", "fa", "txt"])
                separate_file2 = st.file_uploader("Second Sequence FASTA", 
                                               type=["fasta", "fa", "txt"])
            
            if st.button("🧹 Clear Uploads"):
                st.session_state.uploaded_files = None
                st.rerun()
            
            def parse_fasta(uploaded_file):
                """Parse FASTA file and return sequence"""
                if uploaded_file is None:
                    return ""
                    
                try:
                    if uploaded_file.size > 5 * 1024 * 1024:  # 5MB limit
                        st.error("File too large (max 5MB)")
                        return ""
                        
                    file_content = uploaded_file.getvalue().decode("utf-8")
                    records = list(SeqIO.parse(StringIO(file_content), "fasta"))
                    
                    if not records:
                        st.error("No valid sequences found in file")
                        return ""
                        
                    return str(records[0].seq)
                    
                except Exception as e:
                    st.error(f"Error parsing FASTA: {str(e)}")
                    return ""
            
            # Handle combined file upload
            if combined_file:
                try:
                    file_content = combined_file.getvalue().decode("utf-8")
                    records = list(SeqIO.parse(StringIO(file_content), "fasta"))
                    
                    if len(records) >= 2:
                        seq1 = str(records[0].seq)
                        seq2 = str(records[1].seq)
                        
                        with st.expander("👀 Preview Uploaded Sequences"):
                            st.markdown(f"**Sequence 1:** {records[0].description}")
                            st.code(seq1[:100] + ("..." if len(seq1) > 100 else ""))
                            
                            st.markdown(f"**Sequence 2:** {records[1].description}")
                            st.code(seq2[:100] + ("..." if len(seq2) > 100 else ""))
                    else:
                        st.warning("Uploaded file must contain at least 2 sequences")
                        
                except Exception as e:
                    st.error(f"Error processing combined file: {str(e)}")
            
            # Handle separate file uploads
            if separate_file1 and separate_file2:
                seq1 = parse_fasta(separate_file1)
                seq2 = parse_fasta(separate_file2)
                
                if seq1 and seq2:
                    with st.expander("👀 Preview Uploaded Sequences"):
                        st.markdown("**Sequence 1 (First 100 chars):**")
                        st.code(seq1[:100] + ("..." if len(seq1) > 100 else ""))
                        
                        st.markdown("**Sequence 2 (First 100 chars):**")
                        st.code(seq2[:100] + ("..." if len(seq2) > 100 else ""))

        # Alignment parameters
        st.markdown("---")
        st.markdown("### ⚙️ Alignment Parameters")
        col1, col2 = st.columns(2)
        with col1:
            match_score = st.number_input("Match Score", value=1)
            mismatch_score = st.number_input("Mismatch Score", value=-1)
        with col2:
            gap_open = st.number_input("Gap Opening Penalty", value=-0.5)
            gap_extend = st.number_input("Gap Extension Penalty", value=-0.1)

        if st.button("🔗 Perform Alignment", type="primary"):
            if not seq1 or not seq2:
                st.warning("Please provide both sequences before alignment")
            else:
                # Clean sequences (remove whitespace and numbers)
                clean_seq1 = re.sub(r'[^a-zA-Z]', '', seq1).upper()
                clean_seq2 = re.sub(r'[^a-zA-Z]', '', seq2).upper()
                
                # Validate sequences based on type
                valid_chars = set("ACGTU" if seq_type == "Nucleotide" else "ACDEFGHIKLMNPQRSTVWY")
                invalid_chars1 = set(clean_seq1) - valid_chars
                invalid_chars2 = set(clean_seq2) - valid_chars
                
                if invalid_chars1 or invalid_chars2:
                    st.error(f"Invalid characters detected for {seq_type} sequence: "
                            f"{', '.join(invalid_chars1 | invalid_chars2)}")
                else:
                    with st.spinner("Performing alignment..."):
                        try:
                            if alignment_type.startswith("Global"):
                                alignments = pairwise2.align.globalms(
                                    clean_seq1, clean_seq2, 
                                    match_score, mismatch_score, 
                                    gap_open, gap_extend
                                )
                            else:
                                alignments = pairwise2.align.localms(
                                    clean_seq1, clean_seq2,
                                    match_score, mismatch_score,
                                    gap_open, gap_extend
                                )

                            if alignments:
                                alignment_text = ""
                                for i, alignment in enumerate(alignments[:1]):  # Show only top alignment
                                    st.markdown(f"### 🧬 Alignment Result #{i+1}")
                                    formatted_alignment = format_alignment(*alignment)
                                    alignment_text += formatted_alignment + "\n\n"

                                    with st.expander("🔍 View Full Alignment"):
                                        st.code(formatted_alignment, language="text")

                                    # Calculate metrics
                                    score = alignment[2]
                                    aligned_seq1, aligned_seq2 = alignment[0], alignment[1]
                                    matches = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2) if a != '-' and b != '-')
                                    identity = matches / len(aligned_seq1) * 100
                                    gaps = aligned_seq1.count('-') + aligned_seq2.count('-')

                                    # Display metrics
                                    cols = st.columns(4)
                                    cols[0].metric("Alignment Score", f"{score:.2f}")
                                    cols[1].metric("Identity", f"{identity:.1f}%")
                                    cols[2].metric("Alignment Length", len(aligned_seq1))
                                    cols[3].metric("Total Gaps", gaps)

                                # Download option
                                st.download_button(
                                    label="📥 Download Alignment",
                                    data=alignment_text,
                                    file_name="alignment_result.txt",
                                    mime="text/plain"
                                )
                            else:
                                st.warning("No alignment could be generated")
                                
                        except Exception as e:
                            st.error(f"Alignment failed: {str(e)}")

# ------------------------------ SEQUENCE TRANSLATION ------------------------------
elif page == "🧬 Sequence Translation":
    with stylable_container(key="translation_container", css_styles="""
        { background: rgba(255,255,255,0.93); border-radius: 15px; padding: 2rem; }
    """):
        st.title("🧬 DNA → Protein Sequence Translation")
        st.markdown("""
        Translate DNA/RNA sequences to protein sequences using standard genetic codes.
        Supports all NCBI translation tables (bacterial, mitochondrial, etc.).
        """)

        with st.expander("💡 How to use this tool - Click for Examples", expanded=False):
            st.markdown("""
            **Example 1: Basic Translation**
            ```
            Input DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
            Genetic Code: 1 (Standard)
            Output: MAIVMGR*KGAR*
            ```
            
            **Example 2: All Reading Frames**  
            For the same sequence, you'll see:
            ```
            Frame +1: MAIVMGR*KGAR*
            Frame +2: WPL*WGREKVRD
            Frame +3: GHCNGP*KGCPD
            And their reverse complements (-1, -2, -3)
            ```
            
            **Example 3: Mitochondrial Code (Table 2)**  
            ```
            Input DNA: ATGACGTAA
            Table 2 (Vertebrate Mitochondrial)
            Output: MT*
            ```
            """)

        # Sequence input
        dna_seq = st.text_area("Enter DNA/RNA Sequence", height=150,
                             value=st.session_state.get("dna_example", ""),
                             placeholder="ATGCGTACGTAA...",
                             help="Paste your nucleotide sequence (DNA or RNA)")
        
        # Options
        col1, col2 = st.columns(2)
        with col1:
            to_protein = st.radio("Convert to:", ["Protein Sequence", "All Reading Frames"],
                                help="Standard translation or all 6 reading frames")
        with col2:
            table_id = st.selectbox("Genetic Code Table", 
                                  sorted(GENETIC_CODE_NAMES.items()),
                                  format_func=lambda x: f"{x[0]} - {x[1]}",
                                  help="NCBI translation tables")
        
        # Action buttons
        col1, col2 = st.columns([1, 3])
        with col1:
            if st.button("🧬 Translate"):
                if not dna_seq:
                    st.warning("Please enter a DNA/RNA sequence")
                else:
                    # Clean sequence
                    clean_seq = re.sub(r'[^a-zA-Z]', '', dna_seq).upper()
                    
                    try:
                        table_num = table_id[0]
                        table_name = GENETIC_CODE_NAMES.get(table_num, f"Table {table_num}")
                        st.success(f"Using genetic code table {table_num}: {table_name}")
                        
                        # Create sequence object
                        seq_obj = Seq(clean_seq)
                        
                        if to_protein == "Protein Sequence":
                            # Standard translation
                            protein_seq = seq_obj.translate(table=table_num, to_stop=False)
                            
                            st.markdown("### Translation Result")
                            st.code(str(protein_seq), language="text")
                            
                            # Show stats
                            st.markdown("**Sequence Info:**")
                            col1, col2 = st.columns(2)
                            col1.metric("Input Length", f"{len(clean_seq)} nt")
                            col2.metric("Protein Length", f"{len(protein_seq)} aa")
                            
                            # Download button
                            st.markdown(
                                download_link(str(protein_seq), "translated_sequence.faa", 
                                            "📥 Download Protein Sequence (FASTA)"), 
                                unsafe_allow_html=True
                            )
                            
                        else:
                            # All reading frames
                            st.markdown("### All 6 Reading Frames")
                            frames = []
                            
                            for frame in range(3):
                                st.markdown(f"#### Frame +{frame+1}")
                                protein = seq_obj[frame:].translate(table=table_num)
                                st.code(protein, language="text")
                                frames.append(f">Frame+{frame+1}\n{protein}\n")
                                
                                st.markdown(f"#### Frame -{frame+1}")
                                protein = seq_obj.reverse_complement()[frame:].translate(table=table_num)
                                st.code(protein, language="text")
                                frames.append(f">Frame-{frame+1}\n{protein}\n")
                            
                            # Download all frames
                            st.markdown(
                                download_link("\n".join(frames), "all_frames_translation.faa", 
                                            "📥 Download All Frames (FASTA)"), 
                                unsafe_allow_html=True
                            )
                            
                    except Exception as e:
                        st.error(f"Translation error: {str(e)}")
        with col2:
            if st.button("🔄 Load Example Sequence"):
                example_dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
                st.session_state.dna_example = example_dna
                st.rerun()

# ------------------------------ ABOUT ------------------------------
elif page == "📌 About":
    with stylable_container(key="about_container", css_styles="""
        { background: rgba(255,255,255,0.93); border-radius: 15px; padding: 2rem; }
    """):
        st.title("📌 About GeneVista & The Minds Behind It")

        with stylable_container("about-section", css_styles="background-color: #f0f8ff; padding: 1.5rem; border-radius: 10px;"):
            st.image(
                "https://media.licdn.com/dms/image/v2/D4E03AQFIOiRUAc60bQ/profile-displayphoto-shrink_400_400/profile-displayphoto-shrink_400_400/0/1729440388857?e=1749686400&v=beta&t=_cWkoM-scb2h-WwNfpUHliebtwIe1HVnra3uA1okzRE",
                width=150,
                caption="Rushabh Wakade",
                use_container_width=False
            )

            st.markdown("""
            ### 👨‍💻 Author  
            The author is currently pursuing a **Master's degree in Bioinformatics** at **DES Pune University, Maharashtra, India**.  
            This web server, **GeneVista**, was initially developed as a **mini-project** as part of the academic curriculum, with the intention of providing a hands-on, user-friendly tool for sequence analysis in the bioinformatics community.

            The idea behind building this platform was to create an accessible and browser-based solution for essential sequence retrieval, alignment, and analysis tasks without requiring complex installations or advanced programming skills.

            The author holds a strong interest in the areas of **structural bioinformatics**, **genomic data analysis**, and **computational biology tools development**.  
            This project marks an initial step toward developing more advanced, feature-rich bioinformatics applications in the future.

            In addition to this, the author actively participates in academic seminars, workshops, and poster competitions, and has contributed to research initiatives within the field of bioinformatics.

            The web server will continue to be enhanced with additional tools, improved UI, and expanded sequence analysis features based on user feedback and technological trends.

            ---

            ### 🎯 Purpose  
            - To provide an easy-to-use platform for essential bioinformatics tasks.  
            - To simplify sequence retrieval and alignment for students, researchers, and educators.  
            - To demonstrate core biological data handling using Python and bioinformatics libraries.

            ### 🌟 Key Features  
            - 🧬 **NCBI Sequence Fetching** (nucleotide & protein)  
            - 📑 **GenBank Record Viewer** with detailed metadata  
            - 📊 **Annotated Feature Table** with sequence stats  
            - 🔗 **Pairwise Sequence Alignment** (Global & Local)  
            - 🧬 **Sequence Translation** (Standard & All Frames)  
            - ⚙️ **Custom Scoring Parameters** for alignments  
            - 📂 **Downloadable Outputs** (GenBank, TXT alignments, FASTA translations)

            ### 🧰 Tools & Technologies Used  
            - **Python** and **Streamlit** for building the web interface  
            - **Biopython** for NCBI access and sequence operations  
            - **Pandas** for tabular feature rendering  
            - **Streamlit Extras** for UI styling and enhanced interactivity  
            - Standard Python libraries like `base64`, `io`, and `re`

            ### 💡 Benefits  
            - Browser-based access — no installation required  
            - Beginner-friendly design with expandable sections  
            - Fast, integrated access to public biological databases  
            - Modular structure — ready for future upgrades

            ### 🚀 Planned Future Enhancements  
            - Integration of **BLAST search**  
            - Support for **FASTA file uploads**  
            - **GC-content plots** and sequence visualizations  
            - Session-based project saving and report export

            ### Additional Tip 
            - Keep your device settings on Light mode/theme for better visualization of the webserver.

            ---

            ### 👨‍🏫 Mentorship & Acknowledgements  
            Special thanks to **Dr. Kushagra Kashyap**,  
            *Assistant Professor (Bioinformatics)*, Department of Life Sciences,  
            School of Science and Mathematics, **DES Pune University**, for his invaluable guidance and academic support during the development of this project.  

            His mentorship played a crucial role in shaping the project’s scientific direction and refining the technical implementation. The encouragement and expert insights provided throughout the mini-project were instrumental in its successful completion.  

            🔗 [Connect on LinkedIn](https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb)

            ---

            ### ✉️ Feedback & Contact  
            Your feedback is highly appreciated!  
            📧 Email: [wakaderushabh659@gmail.com](mailto:wakaderushabh659@gmail.com)  
            🔗 LinkedIn: [Rushabh Wakade](https://www.linkedin.com/in/rushabh-wakade-624304318)  
            💻 GitHub: [View Source Code](https://github.com/Rushabhbw/GeneVista/blob/main/Web1.py)

            <center>Built with ❤️ by a passionate bioinformatics learner.</center>
            """, unsafe_allow_html=True)

# Footer  
st.markdown("""
<div class="custom-footer">
    <center>2025 Rushabh Wakade | Crafted for Bioinformatics Enthusiasts</center>
</div>
""", unsafe_allow_html=True)
