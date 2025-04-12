import streamlit as st
from Bio import Entrez, SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from io import StringIO
import pandas as pd
import base64
from streamlit_extras.stylable_container import stylable_container
from PIL import Image
import requests
import time
import re

# Set your email for NCBI Entrez
Entrez.email = "wakaderushabh@gmail.com"

# Streamlit page configuration
st.set_page_config(page_title="🔬 GeneVista: Bioinformatics Sequence Hub", layout="wide")

# ---------- Custom CSS for Visual Styling ----------
st.markdown(""" 
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Poppins:wght@300;400;600&display=swap');

    html, body, [class*="css"] {
        font-family: 'Poppins', sans-serif;
    }

    .stApp {
        background: linear-gradient(to right, #eef2f3, #8e9eab);
    }

    .custom-footer {
        margin-top: 50px;
        padding: 10px;
        color: white;
        background-color: #2c3e50;
        border-radius: 10px;
    }

    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
        padding-left: 2rem;
        padding-right: 2rem;
    }

    /* Horizontal Tabs Styling */
    .tab {
        display: flex;
        justify-content: space-around;
        margin-bottom: 20px;
    }

    .tab button {
        background-color: #f1f1f1;
        border: none;
        outline: none;
        padding: 10px 20px;
        cursor: pointer;
        font-size: 18px;
        transition: 0.3s;
        border-radius: 5px 5px 0 0;
        margin: 0 5px; /* Reduced margin to bring tabs closer */
    }

    .tab button:hover {
        background-color: #ddd;
    }

    .tab button.active {
        background-color: #4CAF50;
        color: white;
    }

    /* Animation for Tab Content */
    .tab-content {
        padding: 20px;
        border: 1px solid #ccc;
        border-top: none;
        border-radius: 0 0 5px 5px;
        animation: fadeIn 1s ease-out;
        opacity: 0;
        animation-fill-mode: forwards;
    }

    /* Fade In Effect */
    @keyframes fadeIn {
        0% {
            opacity: 0;
            transform: translateY(10px);
        }
        100% {
            opacity: 1;
            transform: translateY(0);
        }
    }

    /* Tab content visibility */
    .tab-content.show {
        animation: fadeIn 1s ease-out forwards;
    }

    </style>
""", unsafe_allow_html=True)


# JavaScript to force content to show animation on tab click (if necessary):
st.markdown("""
    <script>
    const tabs = document.querySelectorAll('.tab button');
    const contents = document.querySelectorAll('.tab-content');
    
    tabs.forEach((tab, index) => {
        tab.addEventListener('click', () => {
            contents.forEach(content => {
                content.classList.remove('show');
            });
            contents[index].classList.add('show');
        });
    });
    </script>
""", unsafe_allow_html=True)


# --- Top Navigation State Handling ---
if "top_nav" not in st.session_state:
    st.session_state.top_nav = "🏠 Home"

# Render Top Navigation as Buttons
with st.container():
    col1, col2, col3, col4, col5, col6 = st.columns([2, 1, 1, 1, 1, 1])
    col1.markdown("<h1 style='margin-bottom:0;color:#2c3e50;'>GeneVista</h1>", unsafe_allow_html=True)
    if col2.button("🏠 Home"):
        st.session_state.top_nav = "🏠 Home"
    if col3.button("📖 User Guide"):
        st.session_state.top_nav = "📖 User Guide"
    if col4.button("🧬 Fetch Sequence"):
        st.session_state.top_nav = "🧬 Fetch Sequence (NCBI)"
    if col5.button("🔗 Alignment"):
        st.session_state.top_nav = "🔗 Sequence Alignment"
    if col6.button("📌 About"):
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

    - **📑 Downloadable Reports:**  
    Export alignment results and GenBank files for offline analysis and documentation.

    - **📖 Easy-to-Follow User Guide:**  
    A step-by-step manual to help you navigate and use the platform efficiently.

    ---

    ### 🎯 Why Use GeneVista?

    - **No installation needed** — fully browser-based.
    - Beginner-friendly, with expandable sections and interactive controls.
    - Direct access to public biological databases.
    - Designed for learning, research, and practical applications in bioinformatics.

    ---

    💡 *Explore the tools via the sidebar navigation and start analyzing your sequences now!*
    """)

# ------------------------------ USER GUIDE ------------------------------
elif page == "📖 User Guide":
    st.title("📖 User Guide")
    st.markdown("""
    ### 📌 Features Overview:

    1. **🧬 Fetch Sequence (NCBI)**
       - Input either nucleotide or protein accession number.
       - Choose database: `nucleotide` or `protein`.
       - View metadata, description, organism, sequence, and GenBank data.

    2. **🔗 Sequence Alignment**
       - Paste two sequences.
       - Choose between Global (Needleman-Wunsch) and Local (Smith-Waterman) alignment.
       - Choose sequence type: Nucleotide or Protein.
       - Customize scoring and gap penalties.
       - View aligned results including score, identity, and match statistics.

    ### 🔗 Additional Resources:
    - [Biopython Documentation](https://biopython.org/wiki/Documentation)
    - [NCBI Handbook](https://www.ncbi.nlm.nih.gov/books/NBK143764/)
    """)

# ------------------------------ FETCH SEQUENCE ------------------------------
elif page == "🧬 Fetch Sequence (NCBI)":
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
    st.title("🔗 Sequence Alignment")
    alignment_type = st.radio("Choose Alignment Type", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"])
    seq_type = st.radio("Select Sequence Type", ["Nucleotide", "Protein"])

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

    match_score = st.number_input("Match Score", value=1)
    mismatch_score = st.number_input("Mismatch Score", value=-1)
    gap_open = st.number_input("Gap Opening Penalty", value=-0.5)
    gap_extend = st.number_input("Gap Extension Penalty", value=-0.1)

    if st.button("🔗 Perform Alignment"):
        if not seq1 or not seq2:
            st.warning("Please enter both sequences.")
        else:
            if alignment_type.startswith("Global"):
                alignments = pairwise2.align.globalms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)
            else:
                alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)

            if alignments:
                alignment_text = ""
                for i, alignment in enumerate(alignments[:1]):
                    st.markdown(f"### 🧬 Alignment Result #{i+1}")
                    formatted_alignment = format_alignment(*alignment)
                    alignment_text += formatted_alignment + "\n\n"

                    with st.expander("🔍 View Alignment"):
                        st.code(formatted_alignment, language="text")

                    score = alignment[2]
                    start = alignment[3]
                    end = alignment[4]
                    alignment_len = end - start
                    matches = sum(a == b for a, b in zip(alignment[0], alignment[1]) if a != '-' and b != '-')
                    identity_pct = (matches / alignment_len) * 100 if alignment_len > 0 else 0

                    col1, col2, col3, col4 = st.columns(4)
                    col1.metric("Alignment Score", f"{score}")
                    col2.metric("Alignment Length", f"{alignment_len}")
                    col3.metric("Matches", f"{matches}")
                    col4.metric("Identity (%)", f"{identity_pct:.2f}")

                st.markdown("### 📥 Download Alignment Result")
                st.markdown(download_link(alignment_text, "alignment_result.txt", "📂 Download Alignment as TXT"), unsafe_allow_html=True)
            else:
                st.warning("No alignment could be generated.")

# ------------------------------ ABOUT ------------------------------
elif page == "📌 About":
    st.title("📌 About This Bioinformatics Web Server")

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
        - ⚙️ **Custom Scoring Parameters** for alignments  
        - 📂 **Downloadable Outputs** (GenBank, TXT alignments)

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

        ---

        ### ✉️ Feedback & Contact  
        Your feedback is highly appreciated!  
        📧 Email: [wakaderushabh659@gmail.com](mailto:wakaderushabh659@gmail.com)  
        🔗 LinkedIn: [Rushabh Wakade](https://www.linkedin.com/in/rushabh-wakade-624304318)

        <center>Built with ❤️ by a passionate bioinformatics learner.</center>
        """, unsafe_allow_html=True)

# Footer  
st.markdown('<div class="custom-footer"><center>© 2025 Rushabh Wakade | Crafted for Bioinformatics Enthusiasts</center></div>', unsafe_allow_html=True)
