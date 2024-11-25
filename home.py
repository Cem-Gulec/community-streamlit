import streamlit as st

def show_home_page():
    st.markdown("""
    <style>
    .main-title {
        font-size: 36px;
        font-weight: bold;
        margin-bottom: 20px;
    }
    .description {
        font-size: 16px;
        line-height: 1.5;
        margin-bottom: 20px;
    }
    .highlight {
        color: #FF4B4B;
        font-weight: bold;
    }
    .button-container {
        display: flex;
        justify-content: space-between;
        margin-top: 30px;
    }
    .custom-button {
        background-color: #0E1117;
        border: 2px solid white;
        color: #31333F;
        padding: 10px 20px;
        text-align: center;
        text-decoration: none;
        display: inline-block;
        font-size: 16px;
        margin: 4px 2px;
        cursor: pointer;
        border-radius: 5px;
    }
    /* Custom CSS for Streamlit button */
    .stButton > button {
        background-color: #0E1117;
        border: 2px solid white;
        color: #259bc4;
        padding: 10px 20px;
        text-align: center;
        text-decoration: none;
        display: inline-block;
        font-size: 16px;
        margin: 4px 2px;
        cursor: pointer;
        border-radius: 5px;
        width: auto;
        height: auto;
    }
    </style>
    """, unsafe_allow_html=True)

    st.markdown('<div class="main-title">Welcome to Community Web</div>', unsafe_allow_html=True)

    st.markdown("""
    <div class="description">
    <span class="highlight">Community</span> is an R package designed to explore the differences in communication between various case and control samples using single-cell RNA sequencing (scRNAseq). With its user-friendly output,
    <span class="highlight">Community</span> lets you easily follow the overall shifts in communication within the cohorts. You can visualize and delve into the most significant differences in interactions, and even investigate
    what's driving these changes. It's a handy tool for anyone interested in a deeper understanding of cell-to-cell communication.
    </div>
    """, unsafe_allow_html=True)

    col1, col2, col3 = st.columns([4, 5, 5.5])

    with col1:
        st.markdown("""
        <a href="https://github.com/SoloveyMaria/community" target="_blank" class="custom-button">📘 Repository Page</a>
        """, unsafe_allow_html=True)

    with col2:
        if st.button("🧬 Calculate Communication"):
            st.session_state.page = "Perform Analysis"
            st.rerun()

    with col3:
        st.markdown("""
        <a href="https://github.com/SoloveyMaria/community/blob/main/docs/showcase_notebooks/Lasry/calculate_communication.ipynb" target="_blank" class="custom-button">📄 More Documentation</a>
        """, unsafe_allow_html=True)