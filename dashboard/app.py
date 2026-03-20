# dashboard/app.py
"""
pH-ZinCloud — Interactive Zinc-Binding Site Stability Dashboard
Streamlit web application with pH slider and 3D molecular viewer.

Run locally:
    streamlit run dashboard/app.py

Deploy:
    streamlit run dashboard/app.py --server.port 8080
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import requests

# ── Page configuration ─────────────────────────────────────────────────────────
st.set_page_config(
    page_title="pH-ZinCloud",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── Custom CSS ─────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1F4E79;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 0.95rem;
        color: #555;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background: #f0f4f8;
        border-radius: 8px;
        padding: 0.8rem 1rem;
        margin: 0.3rem 0;
    }
    .stable   { color: #2e7d32; font-weight: 600; }
    .moderate { color: #f57c00; font-weight: 600; }
    .unstable { color: #c62828; font-weight: 600; }
    .switch-badge {
        background: #ffebee;
        color: #c62828;
        border-radius: 4px;
        padding: 2px 8px;
        font-size: 0.8rem;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)


# ── Helper functions ───────────────────────────────────────────────────────────

def score_colour(score: float) -> str:
    if score >= 0.75:
        return "#2e7d32"   # green
    elif score >= 0.40:
        return "#f57c00"   # amber
    else:
        return "#c62828"   # red


def score_label(score: float) -> str:
    if score >= 0.75:
        return "High stability"
    elif score >= 0.40:
        return "Moderate stability"
    elif score >= 0.10:
        return "Low stability"
    else:
        return "Disrupted"


def get_ph_score(row: pd.Series, pH: float) -> float:
    """Extract stability score at a given pH from a result row."""
    col = f"pH_{str(pH).replace('.', '_')}_score"
    if col in row.index:
        return float(row[col])
    # Interpolate between nearest available pH values if exact not found
    ph_cols = {
        float(c.replace("pH_", "").replace("_score", "").replace("_", ".")):c
        for c in row.index if c.startswith("pH_") and c.endswith("_score")
    }
    if not ph_cols:
        return 0.0
    nearest = min(ph_cols.keys(), key=lambda x: abs(x - pH))
    return float(row[ph_cols[nearest]])


def build_titration_curve(site_df: pd.DataFrame, zinc_site_id: str) -> go.Figure:
    """Build a Plotly titration curve for one zinc site."""
    site_row = site_df[site_df["zinc_site_id"] == zinc_site_id].iloc[0]

    ph_cols = {
        float(c.replace("pH_", "").replace("_score", "").replace("_", ".")):c
        for c in site_row.index if c.startswith("pH_") and c.endswith("_score")
    }
    ph_values = sorted(ph_cols.keys())
    scores    = [float(site_row[ph_cols[ph]]) for ph in ph_values]

    colours = [score_colour(s) for s in scores]

    fig = go.Figure()

    # Coloured line segments
    for i in range(len(ph_values) - 1):
        fig.add_trace(go.Scatter(
            x=[ph_values[i], ph_values[i+1]],
            y=[scores[i], scores[i+1]],
            mode="lines",
            line=dict(color=colours[i], width=3),
            showlegend=False,
            hoverinfo="skip"
        ))

    # Scatter points
    fig.add_trace(go.Scatter(
        x=ph_values, y=scores,
        mode="markers",
        marker=dict(size=8, color=colours, line=dict(color="white", width=1)),
        text=[f"pH {ph:.1f}: {s:.4f} ({score_label(s)})"
              for ph, s in zip(ph_values, scores)],
        hovertemplate="%{text}<extra></extra>",
        showlegend=False
    ))

    # Threshold lines
    fig.add_hline(y=0.75, line_dash="dot", line_color="#2e7d32",
                  annotation_text="High stability (0.75)",
                  annotation_position="right")
    fig.add_hline(y=0.40, line_dash="dot", line_color="#f57c00",
                  annotation_text="Moderate (0.40)",
                  annotation_position="right")
    fig.add_hline(y=0.10, line_dash="dot", line_color="#c62828",
                  annotation_text="Disrupted (0.10)",
                  annotation_position="right")

    fig.update_layout(
        title=f"Stability Titration Curve — {zinc_site_id}",
        xaxis_title="pH",
        yaxis_title="Site Stability Score",
        yaxis=dict(range=[0, 1.05]),
        xaxis=dict(range=[3.8, 9.2]),
        height=380,
        margin=dict(l=40, r=120, t=50, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    fig.update_xaxes(showgrid=True, gridcolor="#eee")
    fig.update_yaxes(showgrid=True, gridcolor="#eee")

    return fig


def render_3d_viewer(pdb_id: str, site_df: pd.DataFrame,
                     zinc_site_id: str, selected_pH: float):
    """Render 3D molecular viewer using py3Dmol via Streamlit components."""
    import streamlit.components.v1 as components

    try:
        # Fetch PDB structure
        pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        response = requests.get(pdb_url, timeout=15)
        if response.status_code != 200:
            st.warning(f"Could not fetch structure for {pdb_id} from RCSB.")
            return
        pdb_data = response.text

        # Get coordinating residues for this site
        site_residues = site_df[site_df["zinc_site_id"] == zinc_site_id].copy()

        # Build residue colour selections as JavaScript
        residue_styles = []
        for _, res_row in site_residues.iterrows():
            score  = get_ph_score(res_row, selected_pH)
            colour = score_colour(score)
            chain  = res_row["chain"]
            resnum = int(res_row["residue_seq"])
            residue_styles.append(
                f"viewer.addStyle({{chain:'{chain}',resi:{resnum}}},"
                f"{{stick:{{color:'{colour}',radius:0.3}},"
                f"sphere:{{color:'{colour}',radius:0.35}}}});"
            )

        residue_js = "\n".join(residue_styles)

        # Build complete HTML with embedded py3Dmol viewer
        html = f"""
        <html>
        <head>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        </head>
        <body style="margin:0;padding:0;background:#1E2530;">
        <div id="viewer" style="width:100%;height:380px;position:relative;"></div>
        <script>
            $(document).ready(function() {{
                let viewer = $3Dmol.createViewer(
                    document.getElementById('viewer'),
                    {{backgroundColor:'#1E2530'}}
                );

                let pdbData = `{pdb_data.replace('`', "'").replace('\\', '\\\\')}`;

                viewer.addModel(pdbData, 'pdb');
                viewer.setStyle({{}}, {{cartoon:{{color:'#888888',opacity:0.6}}}});

                {residue_js}

                // Zinc as blue sphere
                viewer.addStyle({{resn:'ZN'}},
                    {{sphere:{{color:'#1565c0',radius:0.9}}}});

                viewer.zoomTo({{resn:'ZN'}});
                viewer.render();
            }});
        </script>
        </body>
        </html>
        """

        components.html(html, height=385, scrolling=False)

    except Exception as e:
        st.warning(f"3D viewer error: {e}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN APP
# ══════════════════════════════════════════════════════════════════════════════

def main():
    # ── Header ─────────────────────────────────────────────────────────────────
    st.markdown('<div class="main-header">🧪 pH-ZinCloud</div>',
                unsafe_allow_html=True)
    st.markdown(
        '<div class="sub-header">Cloud-native pH-dependent zinc-binding site '
        'stability analyser — MSc Big Data Analytics, RGU</div>',
        unsafe_allow_html=True
    )

    # ── Sidebar ─────────────────────────────────────────────────────────────────
    with st.sidebar:
        st.header("⚙️ Controls")

        # PDB input
        st.subheader("Protein Input")
        pdb_input = st.text_input(
            "Enter PDB ID",
            value="1CA2",
            max_chars=4,
            help="4-character PDB identifier (e.g. 1CA2, 1ZNF, 3CPA)"
        ).strip().upper()

        analyse_btn = st.button("🔬 Analyse", use_container_width=True, type="primary")

        st.divider()

        # pH slider
        st.subheader("pH Selection")
        selected_pH = st.slider(
            "Select pH",
            min_value=4.0,
            max_value=9.0,
            value=7.4,
            step=0.5,
            help="Adjust pH to see how stability scores change"
        )

        # pH context label
        ph_contexts = {
            4.0: "🔴 Lysosome (pH ~4.5)",
            4.5: "🔴 Lysosome (pH ~4.5)",
            5.0: "🔴 Late endosome (pH ~5.0)",
            5.5: "🟠 Early endosome (pH ~5.5)",
            6.0: "🟠 Tumour microenv. (pH ~6.0)",
            6.5: "🟠 Endosome / tumour (pH ~6.5)",
            7.0: "🟡 Cytosol (pH ~7.2)",
            7.4: "🟢 Blood / physiological (pH 7.4)",
            8.0: "🟢 Slightly alkaline",
            8.5: "🟢 Alkaline",
            9.0: "🟢 Alkaline"
        }
        st.caption(ph_contexts.get(selected_pH, f"pH {selected_pH}"))

        st.divider()

        # Filter options
        st.subheader("Filters")
        show_switches_only = st.checkbox("pH-switch sites only", value=False)
        min_score = st.slider("Min score at selected pH", 0.0, 1.0, 0.0, 0.05)

        st.divider()
        st.caption("pH-ZinCloud v1.0 | Fathima Farwin | RGU MSc 2026")

    # ── Main content ───────────────────────────────────────────────────────────
    from bq_client import get_protein_data, run_pipeline_for_pdb

    # Initialise session state
    if "df" not in st.session_state:
        st.session_state.df = pd.DataFrame()
    if "current_pdb" not in st.session_state:
        st.session_state.current_pdb = ""

    # Load data when Analyse is clicked
    if analyse_btn and pdb_input:
        with st.spinner(f"Loading data for {pdb_input}..."):
            df = get_protein_data(pdb_input)

            if df.empty:
                st.info(
                    f"**{pdb_input}** not found in database. "
                    f"Running pipeline now — this takes ~30 seconds..."
                )
                with st.spinner("Running pH-ZinCloud pipeline..."):
                    df = run_pipeline_for_pdb(pdb_input)

            if df.empty:
                st.error(
                    f"Could not analyse {pdb_input}. "
                    f"Check the PDB ID is valid at rcsb.org"
                )
            else:
                st.session_state.df = df
                st.session_state.current_pdb = pdb_input
                st.success(
                    f"✅ Loaded {pdb_input}: "
                    f"{df['zinc_site_id'].nunique()} zinc site(s), "
                    f"{len(df)} coordinating residues"
                )

    # Display results if data is loaded
    df = st.session_state.df
    current_pdb = st.session_state.current_pdb

    if df.empty:
        # Welcome screen
        st.info(
            "👈 Enter a PDB ID in the sidebar and click **Analyse** to begin.\n\n"
            "**Try these examples:**\n"
            "- `1CA2` — Carbonic Anhydrase II (3-His catalytic zinc)\n"
            "- `1ZNF` — Classical Cys₂His₂ zinc finger\n"
            "- `3CPA` — Carboxypeptidase A (His₂Glu site)\n"
            "- `4TLN` — Thermolysin (His₂Glu catalytic zinc)"
        )
        return

    # ── Apply filters ──────────────────────────────────────────────────────────
    display_df = df.copy()
    if show_switches_only:
        switch_sites = display_df[display_df["is_ph_switch"] == True]["zinc_site_id"].unique()
        display_df = display_df[display_df["zinc_site_id"].isin(switch_sites)]

    # Get unique sites
    site_ids = display_df["zinc_site_id"].unique().tolist()

    if not site_ids:
        st.warning("No sites match the current filters.")
        return

    # ── Site selector ──────────────────────────────────────────────────────────
    if len(site_ids) > 1:
        selected_site = st.selectbox(
            "Select zinc site to visualise",
            options=site_ids,
            format_func=lambda x: x
        )
    else:
        selected_site = site_ids[0]

    # ── Main layout: 3 columns ─────────────────────────────────────────────────
    col_viewer, col_chart, col_table = st.columns([1.2, 1.1, 0.9])

    # ── Column 1: 3D Viewer ────────────────────────────────────────────────────
    with col_viewer:
        st.subheader("🔬 3D Structure")

        # Stability at current pH
        site_row = display_df[
            display_df["zinc_site_id"] == selected_site
        ].iloc[0]
        current_score = get_ph_score(site_row, selected_pH)
        colour = score_colour(current_score)

        st.markdown(
            f"**{selected_site}** &nbsp;|&nbsp; "
            f"pH {selected_pH} score: "
            f"<span style='color:{colour};font-weight:700'>"
            f"{current_score:.4f} — {score_label(current_score)}"
            f"</span>",
            unsafe_allow_html=True
        )

        if site_row.get("is_ph_switch", False):
            st.markdown('<span class="switch-badge">⚡ pH-Switch Candidate</span>',
                        unsafe_allow_html=True)

        render_3d_viewer(current_pdb, display_df, selected_site, selected_pH)

        # Colour legend
        st.caption(
            "🟢 Stable (≥0.75) &nbsp; 🟠 Moderate (0.40–0.75) &nbsp; "
            "🔴 Disrupted (<0.40)"
        )

    # ── Column 2: Titration curve ──────────────────────────────────────────────
    with col_chart:
        st.subheader("📈 Stability vs pH")
        fig = build_titration_curve(display_df, selected_site)
        st.plotly_chart(fig, use_container_width=True)

        # Residue-level scores at selected pH
        st.subheader("🔑 Coordinating Residues")
        site_residues = display_df[
            display_df["zinc_site_id"] == selected_site
        ][["residue_name", "residue_seq", "chain",
           "pka_value", "pka_tier", "distance_to_zinc"]].copy()

        site_residues["score@pH"] = site_residues.apply(
            lambda r: get_ph_score(
                display_df[
                    (display_df["zinc_site_id"] == selected_site) &
                    (display_df["residue_seq"] == r["residue_seq"])
                ].iloc[0],
                selected_pH
            ), axis=1
        )
        site_residues["pKa source"] = site_residues["pka_tier"].map(
            {1: "Expt.", 2: "PROPKA", 3: "Canonical"}
        )
        site_residues = site_residues.rename(columns={
            "residue_name": "Residue",
            "residue_seq":  "Seq",
            "chain":        "Chain",
            "pka_value":    "pKa",
            "distance_to_zinc": "Dist(Å)"
        })

        st.dataframe(
            site_residues[["Residue","Seq","Chain","pKa","pKa source","Dist(Å)","score@pH"]],
            use_container_width=True,
            hide_index=True
        )

    # ── Column 3: All sites summary ────────────────────────────────────────────
    with col_table:
        st.subheader("📊 All Sites")

        # Build summary table
        summary_rows = []
        for site_id in site_ids:
            site_data = display_df[display_df["zinc_site_id"] == site_id]
            if site_data.empty:
                continue
            row = site_data.iloc[0]
            ph_score = get_ph_score(row, selected_pH)

            # Filter by minimum score
            if ph_score < min_score:
                continue

            summary_rows.append({
                "Site":        site_id.replace(f"{current_pdb}_", ""),
                "Type":        row.get("zinc_site_id", "").split("_")[2] if "_" in str(row.get("zinc_site_id","")) else "—",
                f"Score@{selected_pH}": round(ph_score, 4),
                "Label":       score_label(ph_score),
                "pH-Switch":   "⚡ YES" if row.get("is_ph_switch", False) else "—",
                "Tier":        f"T{row.get('confidence_tier', '?')}",
            })

        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            st.dataframe(summary_df, use_container_width=True, hide_index=True)

            # CSV download
            csv = display_df.to_csv(index=False)
            st.download_button(
                label="⬇️ Download Full Results CSV",
                data=csv,
                file_name=f"{current_pdb}_phzincloud_results.csv",
                mime="text/csv",
                use_container_width=True
            )
        else:
            st.info("No sites match current filters.")

        # Dataset statistics
        st.subheader("📉 Score Distribution")
        ph_score_col = f"pH_{str(selected_pH).replace('.','_')}_score"
        if ph_score_col in display_df.columns:
            site_scores = display_df.groupby("zinc_site_id")[ph_score_col].first()
            fig_hist = go.Figure(go.Histogram(
                x=site_scores.values,
                nbinsx=10,
                marker_color="#2E75B6",
                opacity=0.8
            ))
            fig_hist.update_layout(
                xaxis_title="Stability Score",
                yaxis_title="Count",
                height=200,
                margin=dict(l=20, r=20, t=20, b=30),
                plot_bgcolor="white"
            )
            st.plotly_chart(fig_hist, use_container_width=True)


if __name__ == "__main__":
    main()