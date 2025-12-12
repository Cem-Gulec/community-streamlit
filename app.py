import numpy as np
import os
import pandas as pd
import plotly.graph_objs as go
import plotly.subplots as sp
import streamlit as st
from streamlit_modal import Modal
import subprocess
import time
from PIL import Image
from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
from st_aggrid.shared import JsCode
from st_aggrid import ColumnsAutoSizeMode
from pathlib import Path

from home import show_home_page
from license import show_license_page


def prepare_comm_result(interactions_df):
    weights = interactions_df[['log2FC_weights']].abs()
    
    anno_samples = pd.DataFrame({'case_or_control': ['case', 'control']})
    
    anno_interactions = interactions_df.copy()
    
    return {
        'weights': weights,
        'anno_samples': anno_samples,
        'anno_interactions': anno_interactions
    }


def filter_interactions(comm_result, threshold_log10_cum_weight=0.05, 
                        threshold_frac_samples_per_condition=0.8,
                        threshold_log10_meanexpr_per_condition=0.1, verbose=True):
    comm_result['anno_interactions']['log10_cum_weight'] = np.log10(comm_result['weights']['log2FC_weights'].abs() + 1)
    
    comm_result['anno_interactions']['frac_samples_controls'] = (comm_result['weights']['log2FC_weights'] != 0).astype(float)
    comm_result['anno_interactions']['frac_samples_cases'] = (comm_result['weights']['log2FC_weights'] != 0).astype(float)

    # Store thresholds
    comm_result['thresholds'] = {
        'threshold_log10_cum_weight': threshold_log10_cum_weight,
        'threshold_frac_samples_per_condition': threshold_frac_samples_per_condition,
        'threshold_log10_meanexpr_per_condition': threshold_log10_meanexpr_per_condition
    }

    # Create plots
    fig_cumW = plot_cumW(comm_result['anno_interactions'], threshold_log10_cum_weight)
    fig_fracSamp = plot_fracSamples(comm_result['anno_interactions'], threshold_frac_samples_per_condition)
    fig_meanLigRec = plot_meanLig_vs_meanRec(comm_result['anno_interactions'], threshold_log10_meanexpr_per_condition)

    # Apply filters
    comm_result['anno_interactions']['passed_log10_cum_weight_filter'] = comm_result['anno_interactions']['log10_cum_weight'] > threshold_log10_cum_weight
    comm_result['anno_interactions']['passed_frac_samples_filter'] = (comm_result['anno_interactions']['frac_samples_controls'] > threshold_frac_samples_per_condition) | \
                                                                     (comm_result['anno_interactions']['frac_samples_cases'] > threshold_frac_samples_per_condition)
    comm_result['anno_interactions']['passed_log10_meanexpr_control_filter'] = (np.log10(comm_result['anno_interactions']['mean_e_s_l_control'] + 1) > threshold_log10_meanexpr_per_condition) & \
                                                                               (np.log10(comm_result['anno_interactions']['mean_e_r_r_control'] + 1) > threshold_log10_meanexpr_per_condition)
    comm_result['anno_interactions']['passed_log10_meanexpr_case_filter'] = (np.log10(comm_result['anno_interactions']['mean_e_s_l_case'] + 1) > threshold_log10_meanexpr_per_condition) & \
                                                                            (np.log10(comm_result['anno_interactions']['mean_e_r_r_case'] + 1) > threshold_log10_meanexpr_per_condition)
    comm_result['anno_interactions']['passed_log10_meanexpr_per_condition_filter'] = comm_result['anno_interactions']['passed_log10_meanexpr_control_filter'] | \
                                                                                     comm_result['anno_interactions']['passed_log10_meanexpr_case_filter']
    comm_result['anno_interactions']['passed_QC_filter'] = (comm_result['anno_interactions']['passed_log10_cum_weight_filter'] & 
                                                            comm_result['anno_interactions']['passed_frac_samples_filter'] & 
                                                            comm_result['anno_interactions']['passed_log10_meanexpr_per_condition_filter'])

    if verbose:
        print(f"{sum(~(comm_result['anno_interactions']['passed_log10_cum_weight_filter'] & comm_result['anno_interactions']['passed_frac_samples_filter']))} out of {len(comm_result['weights'])} interactions do not pass the thresholds for log10 cumulative interactions weight > {threshold_log10_cum_weight} and fraction of expressing samples > {threshold_frac_samples_per_condition}. Also {sum(~comm_result['anno_interactions']['passed_log10_meanexpr_per_condition_filter'])} interactions didn't pass the discrepancy filter. In total, {sum(~comm_result['anno_interactions']['passed_QC_filter'])} bad quality interactions will be removed and {sum(comm_result['anno_interactions']['passed_QC_filter'])} good quality interactions will remain.")

    return comm_result, fig_cumW, fig_fracSamp, fig_meanLigRec


def plot_cumW(df, threshold_log10_cum_weight):
    fig = go.Figure()
    
    # Create histogram trace
    hist_trace = go.Histogram(x=df['log10_cum_weight'], name='Density', yaxis='y2', 
                              marker_color='rgba(0,0,0,0.7)', opacity=0.7)
    
    # Create cumulative distribution trace
    sorted_weights = np.sort(df['log10_cum_weight'])
    cumulative = np.arange(1, len(sorted_weights) + 1) / len(sorted_weights)
    cumulative_trace = go.Scatter(x=sorted_weights, y=cumulative, name='Cumulative', yaxis='y')

    fig.add_trace(hist_trace)
    fig.add_trace(cumulative_trace)

    # Add threshold line
    fig.add_vline(x=threshold_log10_cum_weight, line_dash="dash", line_color="red")

    # Update layout
    fig.update_layout(
        title="Interaction Weight Filter",
        xaxis_title="log10_cum_weight",
        yaxis_title="Cumulative Fraction",
        yaxis2=dict(title="Density", overlaying='y', side='right'),
        showlegend=False
    )

    return fig


def plot_fracSamples(df, threshold_frac_samples_per_condition):
    fig = sp.make_subplots(rows=2, cols=2, 
                           column_widths=[0.2, 0.8], 
                           row_heights=[0.8, 0.2],
                           specs=[[{"type": "histogram"}, {"type": "scatter"}],
                                  [None, {"type": "histogram"}]])

    # Main scatter plot
    fig.add_trace(go.Scatter(
        x=df['frac_samples_controls'],
        y=df['frac_samples_cases'],
        mode='markers',
        marker=dict(
            color='rgba(0,0,0,0.5)',
            size=5
        ),
        name='Interactions'
    ), row=1, col=2)

    # Add threshold lines
    fig.add_hline(y=threshold_frac_samples_per_condition, line_dash="dash", line_color="red", row=1, col=2)
    fig.add_vline(x=threshold_frac_samples_per_condition, line_dash="dash", line_color="red", row=1, col=2)

    # Histogram for x-axis (bottom)
    fig.add_trace(go.Histogram(x=df['frac_samples_controls'], name='Controls', marker_color='rgba(0,0,0,0.5)'), row=2, col=2)

    # Histogram for y-axis (left)
    fig.add_trace(go.Histogram(y=df['frac_samples_cases'], name='Cases', marker_color='rgba(0,0,0,0.5)'), row=1, col=1)

    fig.update_layout(
        title="Fraction of Samples in Which an Interaction is Detected",
        xaxis_title="Controls",
        yaxis_title="Cases",
        showlegend=False
    )

    return fig


def plot_meanLig_vs_meanRec(df, threshold_log10_meanexpr_per_condition):
    fig = sp.make_subplots(rows=1, cols=2, subplot_titles=("Control", "Case"))

    for i, condition in enumerate(['control', 'case']):
        fig.add_trace(go.Scatter(
            x=np.log10(df[f'mean_e_s_l_{condition}'] + 1),
            y=np.log10(df[f'mean_e_r_r_{condition}'] + 1),
            mode='markers',
            name=condition.capitalize()
        ), row=1, col=i+1)

        fig.add_hline(y=threshold_log10_meanexpr_per_condition, line_dash="dash", line_color="red", row=1, col=i+1)
        fig.add_vline(x=threshold_log10_meanexpr_per_condition, line_dash="dash", line_color="red", row=1, col=i+1)

    fig.update_layout(title="Mean Ligand vs Mean Receptor Expression",
                      xaxis_title="log10(mean ligand expression + 1)",
                      yaxis_title="log10(mean receptor expression + 1)")
    return fig


def plot_aml_healthy_visual(aml_data, healthy_data):
    # Create Plotly figure for AML data
    fig_aml = go.Figure()
    fig_aml.add_trace(go.Scatter(
        x=aml_data["number_of_interactions"],
        y=aml_data["mean_interaction_weight"],
        mode="markers",
        marker=dict(
            color=aml_data["interaction_type"].map({
                "among immune cells": "#1f77b4",
                "engages Ery": "#ff7f0e",
                "engages Gran": "#9467bd"
            }),
            size=aml_data["mean_interaction_weight"] * 1000
        ),
        text=aml_data["interaction_ID"]
    ))
    fig_aml.update_layout(
        title="AML",
        xaxis_title="number of interactions",
        yaxis_title="log10 mean w",
        plot_bgcolor="#ffffff",
        paper_bgcolor="#1e1e1e",
        font_color="white"
    )

    # Create Plotly figure for Healthy data
    fig_healthy = go.Figure()
    fig_healthy.add_trace(go.Scatter(
        x=healthy_data["number_of_interactions"],
        y=healthy_data["mean_interaction_weight"],
        mode="markers",
        marker=dict(
            color=healthy_data["interaction_type"].map({
                "among immune cells": "#1f77b4",
                "engages Ery": "#ff7f0e",
                "engages Gran": "#9467bd"
            }),
            size=healthy_data["mean_interaction_weight"] * 1000
        ),
        text=healthy_data["interaction_ID"]
    ))
    fig_healthy.update_layout(
        title="Healthy",
        xaxis_title="number of interactions",
        yaxis_title="log10 mean w",
        plot_bgcolor="#ffffff",
        paper_bgcolor="#1e1e1e",
        font_color="white"
    )

    return fig_aml, fig_healthy


def update_network_plots(immune_color, ery_color, gran_color):
    """Update network plots with new colors by calling R function"""
    # Create a temporary R script
    temp_r_script = """
    source("visualization.R")
    update_network_colors(
        immune_color = "{immune_color}",
        ery_color = "{ery_color}",
        gran_color = "{gran_color}"
    )
    """.format(
        immune_color=immune_color,
        ery_color=ery_color,
        gran_color=gran_color
    )
    
    # Write the temporary script
    with open("temp_update_colors.R", "w") as f:
        f.write(temp_r_script)
    
    try:
        # Run the temporary script
        result = subprocess.run(
            ["Rscript", "temp_update_colors.R"], 
            capture_output=True,
            text=True,
            check=True
        )
        
        if result.stdout:
            print("R Output:", result.stdout)
        if result.stderr:
            print("R Error:", result.stderr)
            
        os.remove("temp_update_colors.R")
        return True
        
    except subprocess.CalledProcessError as e:
        st.error(f"Error updating network plots: {str(e)}")
        print("R Output:", e.stdout)
        print("R Error:", e.stderr)
        # Clean up even if there's an error
        if os.path.exists("temp_update_colors.R"):
            os.remove("temp_update_colors.R")
        return False


def get_top_interactions(n_top=10, interaction_type="all", sort_by="log2FC"):
    """
    Get top n interactions from saved data
    interaction_type: "all", "up", or "down"
    sort_by: "log2FC", "mean_expression", or "p_value"
    """
    try:
        # Load the interactions data with tab separator
        interactions_df = pd.read_csv("visualizations/anno_interactions.txt", sep="\t")
        
        # Calculate absolute values for sorting
        interactions_df['abs_log2FC_weights'] = abs(interactions_df['log2FC_weights'])
        interactions_df['mean_expression_change'] = (
            abs(interactions_df['mean_e_s_l_case'] - interactions_df['mean_e_s_l_control']) +
            abs(interactions_df['mean_e_r_r_case'] - interactions_df['mean_e_r_r_control'])
        ) / 2
        
        if interaction_type == "up":
            filtered_df = interactions_df[interactions_df['log2FC_weights'] > 0]
        elif interaction_type == "down":
            filtered_df = interactions_df[interactions_df['log2FC_weights'] < 0]
        else:
            filtered_df = interactions_df
            
        # Sort based on user selection
        if sort_by == "log2FC":
            sorted_df = filtered_df.sort_values(by=['abs_log2FC_weights'], ascending=[False])
        elif sort_by == "mean_expression":
            sorted_df = filtered_df.sort_values(by=['mean_expression_change'], ascending=[False])
        elif sort_by == "p_value":
            sorted_df = filtered_df.sort_values(by=['p.value'], ascending=[False])
        
        # Select relevant columns
        display_columns = [
            'interaction_ID',
            'ligand_gene_name',
            'receptor_gene_name',
            'sending_cell_type',
            'receiving_cell_type',
            'log2FC_weights',
            'mean_e_s_l_control',
            'mean_e_s_l_case',
            'mean_e_r_r_control',
            'mean_e_r_r_case',
            'mean_expression_change',
            'p.value',
            'p.adj',
            'interaction_category'
        ]
        
        result_df = sorted_df[display_columns].head(n_top)
        
        # Rename columns for better display
        result_df = result_df.rename(columns={
            'ligand_gene_name': 'Ligand',
            'receptor_gene_name': 'Receptor',
            'sending_cell_type': 'Sending Cell',
            'receiving_cell_type': 'Receiving Cell',
            'log2FC_weights': 'Log2FC',
            'mean_e_s_l_control': 'Mean Ligand Expr Control',
            'mean_e_s_l_case': 'Mean Ligand Expr Case',
            'mean_e_r_r_control': 'Mean Receptor Expr Control',
            'mean_e_r_r_case': 'Mean Receptor Expr Case',
            'mean_expression_change': 'Mean Expression Change',
            'p.value': 'P-value',
            'p.adj': 'Adjusted P-value',
            'interaction_category': 'Interaction Category'
        })
        
        return result_df
        
    except Exception as e:
        print(f"Error getting top interactions: {str(e)}")
        st.error(f"Error loading interactions: {str(e)}")
        return pd.DataFrame()


def display_visualizations():
    if 'show_more_forest_plots' not in st.session_state:
        st.session_state.show_more_forest_plots = False
    
    if 'network_colors' not in st.session_state:
        st.session_state.network_colors = {
            'immune': '#1f77b4',  # default blue
            'ery': '#ff7f0e',     # default orange
            'gran': '#9467bd'     # default purple
    }
        
    # Initialize heatmap-specific session states
    if 'selected_interactions' not in st.session_state:
        st.session_state.selected_interactions = [
            "T:CALM1_HSPC:PDE1B", "T:CALM2_HSPC:PDE1B", "HSPC:CALM2_HSPC:PDE1B",
            "HSPC:CALM1_HSPC:PDE1B", "HSPC:CD47_HSPC:SIRPA", "B:CALM2_HSPC:PDE1B",
            "B:CALM1_HSPC:PDE1B", "Mono:CALM1_HSPC:PDE1B", "Mono:CALM3_HSPC:PDE1B",
            "Mono:NRG1_Gran:NETO2", "Mono:NRG1_Gran:MS4A4A", "Mono:NRG1_Gran:HLA-DPB1",
            "Mono:TNFSF14_Gran:LTBR", "Mono:TNFSF14_Gran:TNFRSF14",
            "Mono:TNFSF14_Mono:TNFRSF14", "Mono:NRG1_Mono:SIGLEC7", "Mono:HP_Gran:CD163",
            "Mono:HP_Gran:ASGR2", "Mono:HP_Gran:ASGR1"
        ]
    if 'heatmap_generated' not in st.session_state:
        st.session_state.heatmap_generated = False
    
    if 'selected_interactions' not in st.session_state:
        st.session_state.selected_interactions = []

    image_dir = "./plots/"
    forest_images_dir = "./plots/forest/"
    forest_image_order = [
        "forestplot_no_change.png",
        "forestplot_only_down.png",
        "forestplot_only_up.png",
        "forestplot_concordantDown_one_several.png",
        "forestplot_concordantUp_one_several.png",
        "forestplot_insuffDown.png",
        "forestplot_insuffUp.png",
        "forestplot_suffComp.png"
    ]
   
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Mean weight vs Mean number", 
        "Volcano", 
        "Heatmap", 
        "Network plots", 
        "Forest plots"])

    tab1.subheader("Plotting mumber of interactions vs mean interaction weights per cell type to cell type interaction")
    tab2.subheader("Visualization of differential interactions as a Volcano plot")
    tab3.subheader("Heatmap of interactions weight of top differential interactions")
    tab4.subheader("Visualization of differential interactions as Network plots")
    tab5.subheader("Visualization of individual components")
    
    with tab1:
        col1, col2 = st.columns(2)
        with col1:
            load_and_display_image(os.path.join(image_dir, "Log10_vs_NrInteractions_plot_1.png"))
        with col2:
            load_and_display_image(os.path.join(image_dir, "Log10_vs_NrInteractions_plot_2.png"))
    
    with tab2:
        col1, col2 = st.columns(2)
        
        with col1:
            load_and_display_image(os.path.join(image_dir, "Volcano_plot_.png"))
        
        st.write("### Top Interactions")
        
        # Add controls for the table
        col1, col2, col3 = st.columns(3)
        
        with col1:
            interaction_type = st.radio(
                "Select interaction type:",
                ["All", "Upregulated", "Downregulated"],
                key="volcano_interaction_type"
            )
        
        with col2:
            sort_by = st.radio(
                "Sort by:",
                ["Log2 Fold Change", "Mean Expression Change", "P-value"],
                key="volcano_sort_by"
            )
        
        with col3:
            n_top = st.slider(
                "Number of top interactions to show:",
                min_value=5,
                max_value=50,
                value=10,
                step=5,
                key="volcano_n_top"
            )
        
        # Convert selections to parameters
        interaction_type_param = interaction_type.lower()
        if interaction_type_param == "upregulated":
            interaction_type_param = "up"
        elif interaction_type_param == "downregulated":
            interaction_type_param = "down"
        else:
            interaction_type_param = "all"
            
        if sort_by == "Log2 Fold Change":
            sort_by_param = "log2FC"
        elif sort_by == "Mean Expression Change":
            sort_by_param = "mean_expression"
        else:
            sort_by_param = "p_value"
        
        # Get and display top interactions
        top_df = get_top_interactions(n_top, interaction_type_param, sort_by_param)
        
        if not top_df.empty:
            # Configure grid options
            gb = GridOptionsBuilder.from_dataframe(top_df)
            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
            gb.configure_default_column(
                groupable=True,
                value=True,
                enableRowGroup=True,
                aggFunc="sum",
                editable=False
            )
            
            # Format numeric columns
            numeric_columns = [
                "Log2FC", 
                "Mean Ligand Expr Control", 
                "Mean Ligand Expr Case",
                "Mean Receptor Expr Control",
                "Mean Receptor Expr Case",
                "Mean Expression Change",
                "P-value",
                "Adjusted P-value"
            ]
            
            for col in numeric_columns:
                gb.configure_column(
                    col,
                    type=["numericColumn", "numberColumnFilter"],
                    precision=3
                )
            
            # Configure specific formatting for p-values (scientific notation)
            for col in ["P-value", "Adjusted P-value"]:
                gb.configure_column(
                    col,
                    type=["numericColumn", "numberColumnFilter"],
                    valueFormatter="data < 0.001 ? data.toExponential(3) : data.toFixed(3)"
                )
            
            gridOptions = gb.build()
            
            # Display the grid
            AgGrid(
                top_df,
                gridOptions=gridOptions,
                data_return_mode='AS_INPUT', 
                update_mode=GridUpdateMode.MODEL_CHANGED,
                fit_columns_on_grid_load=True,
                theme='streamlit',
                height=400,
                width='100%',
                reload_data=True,
                columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS
            )
            
            # Add download button for the table
            csv = top_df.to_csv(index=False)
            st.download_button(
                label="Download table as CSV",
                data=csv,
                file_name="top_interactions.csv",
                mime="text/csv"
            )
        else:
            st.warning("No interaction data available")
  
    with tab3:
        try:
            # Read interactions file
            interactions_df = pd.read_csv("visualizations/anno_interactions.txt", sep="\t")
            
            # Create a form to contain the multiselect and generate button
            with st.form("heatmap_form"):
                # Create multiselect widget inside the form
                selected_ids = st.multiselect(
                    "Select Interaction IDs for Heatmap",
                    options=interactions_df['interaction_ID'].tolist(),
                    default=st.session_state.selected_interactions,
                    key="heatmap_interaction_ids"
                )
                
                # Add the generate button as a form submit button
                generate_pressed = st.form_submit_button("Generate Heatmap")
                
                if generate_pressed:
                    if selected_ids:
                        with st.spinner("Generating heatmap..."):
                            # Update session state
                            st.session_state.selected_interactions = selected_ids
                            
                            # Join the selected IDs into a single string
                            string_arg = " ".join(selected_ids)
                            
                            # Run the R script
                            try:
                                result = subprocess.run(
                                    ["Rscript", "heatmap.R", string_arg],
                                    capture_output=True,
                                    text=True,
                                    check=True
                                )
                                
                                if result.stdout:
                                    st.write("R Output:", result.stdout)
                                
                                st.session_state.heatmap_generated = True
                                
                            except subprocess.CalledProcessError as e:
                                st.error(f"Error running R script: {e}")
                                if e.stderr:
                                    st.error(f"R Error Output: {e.stderr}")
                    else:
                        st.warning("Please select at least one interaction ID.")
            
            # Display existing heatmap if it exists
            if os.path.exists(os.path.join("plots", "heatmap.png")):
                if not st.session_state.heatmap_generated:
                    st.info("Select interactions and click 'Generate Heatmap' to update the visualization.")
                
                col1, col2 = st.columns(2)
                with col1:
                    load_and_display_image(os.path.join("plots", "heatmap.png"))
                
        except FileNotFoundError:
            st.error("Interactions file not found. Please ensure anno_interactions.txt is present in the visualizations directory.")
        except Exception as e:
            st.error(f"An error occurred: {e}")


    with tab4:
        col1, col2, col3 = st.columns(3)
            
        with col1:
            new_immune_color = st.color_picker(
                "Immune Cells Color",
                st.session_state.network_colors['immune']
            )
        with col2:
            new_ery_color = st.color_picker(
                "Engages Ery Color",
                st.session_state.network_colors['ery']
            )
        with col3:
            new_gran_color = st.color_picker(
                "Engages Gran Color",
                st.session_state.network_colors['gran']
            )
        
        # Only update if colors have changed
        st.session_state.colors_changed = (
            new_immune_color != st.session_state.network_colors['immune'] or
            new_ery_color != st.session_state.network_colors['ery'] or
            new_gran_color != st.session_state.network_colors['gran']
        )
        
        if st.session_state.colors_changed:
            if st.button("Apply New Colors"):
                with st.spinner("Updating network plots..."):
                    if update_network_plots(new_immune_color, new_ery_color, new_gran_color):
                        # Update session state colors
                        st.session_state.network_colors['immune'] = new_immune_color
                        st.session_state.network_colors['ery'] = new_ery_color
                        st.session_state.network_colors['gran'] = new_gran_color
                        st.success("Network plots updated successfully!")
                        time.sleep(1)
                        st.experimental_rerun()
                    else:
                        st.error("Failed to update network plots. Check the console for details.")
    
    
        # Display network plots
        col1, col2 = st.columns(2)
        with col1:
            load_and_display_image(os.path.join(image_dir, "network_plot_AML.png"))
        with col2:
            load_and_display_image(os.path.join(image_dir, "network_plot_healthy.png"))

    with tab5:
        col1, col2 = st.columns(2)
        with col1:
            if len(forest_image_order) > 0:
                image_name = forest_image_order[0]
                load_and_display_image(os.path.join(forest_images_dir, image_name), 
                                       caption=image_name[11:-4].replace('_', ' ').title())
        with col2:
            if len(forest_image_order) > 1:
                image_name = forest_image_order[1]
                load_and_display_image(os.path.join(forest_images_dir, image_name), 
                                       caption=image_name[11:-4].replace('_', ' ').title())
        
        if len(forest_image_order) > 2:
            if not st.session_state.show_more_forest_plots:
                if st.button("Show More Forest Plots"):
                    st.session_state.show_more_forest_plots = True
                    st.experimental_rerun()
            
            if st.session_state.show_more_forest_plots:
                for i in range(2, len(forest_image_order), 2):
                    col1, col2 = st.columns(2)
                    with col1:
                        if i < len(forest_image_order):
                            image_name = forest_image_order[i]
                            load_and_display_image(os.path.join(forest_images_dir, image_name), 
                                                   caption=image_name[11:-4].replace('_', ' ').title())
                    with col2:
                        if i+1 < len(forest_image_order):
                            image_name = forest_image_order[i+1]
                            load_and_display_image(os.path.join(forest_images_dir, image_name), 
                                                   caption=image_name[11:-4].replace('_', ' ').title())

    
def load_and_display_image(image_path, caption=None):
    try:
        image = Image.open(image_path)
        if caption:
            st.caption(caption)
        st.image(image, use_column_width=True)
    except FileNotFoundError:
        st.error(f"Image not found: {image_path}")
    except Exception as e:
        st.error(f"Error loading image {image_path}: {e}")


def load_data(file):
    # If file is a string (file path)
    if isinstance(file, str):
        file_path = file
        file_name = os.path.basename(file_path)
    
    # If file is a file object (from file uploader)
    else:
        file_path = file
        file_name = file.name

    if file_name.endswith('.csv'):
        df = pd.read_csv(file_path)
    else:
        df = pd.read_csv(file_path, sep='\t', quotechar='"')
    return df


def load_default_files():
    # Get the directory of the current script
    current_dir = Path(__file__).parent.absolute()
    
    # Construct paths to input files
    default_counts_file = current_dir / "input_data" / "toy_counts.csv"
    default_cell_annot_file = current_dir / "input_data" / "toy_cell_annot.csv"
    default_sample_annot_file = current_dir / "input_data" / "anno_samples_norm.txt"
    default_database_file = current_dir / "input_data" / "ready_database.txt"

    # Check if files exist and load them
    if default_counts_file.exists():
        st.session_state.df_counts = load_data(str(default_counts_file))
    else:
        st.error(f"File not found: {default_counts_file}")

    if default_cell_annot_file.exists():
        st.session_state.df_cell_annot = load_data(str(default_cell_annot_file))
    else:
        st.error(f"File not found: {default_cell_annot_file}")

    if default_sample_annot_file.exists():
        st.session_state.df_sample_annot = load_data(str(default_sample_annot_file))
    else:
        st.error(f"File not found: {default_sample_annot_file}")

    if default_database_file.exists():
        st.session_state.df_database = load_data(str(default_database_file))
    else:
        st.error(f"File not found: {default_database_file}")

    if all(file.exists() for file in [default_counts_file, default_cell_annot_file, default_sample_annot_file, default_database_file]):
        st.session_state.files_uploaded = True
    else:
        st.warning("Not all demo dataset files were found. Please check the input_data directory.")


def run_r_script(filename):
    try:
        subprocess.run(["Rscript", filename], check=True)
        return True
    except subprocess.CalledProcessError as e:
        st.error(f"Error running R script: {e}")
        return False


def perform_analysis_callback():
    st.session_state.analysis_performed = True


def perform_analysis():
    st.title("Community")

    # Initializing session states
    if 'files_uploaded' not in st.session_state:
        st.session_state.files_uploaded = False
    if 'df_sample_annot' not in st.session_state:
        st.session_state.df_sample_annot = None
    if 'df_counts' not in st.session_state:
        st.session_state.df_counts = None
    if 'df_database' not in st.session_state:
        st.session_state.df_database = None
    if 'analysis_performed' not in st.session_state: 
        st.session_state.analysis_performed = False
    if 'show_popup' not in st.session_state:
        st.session_state.show_popup = False
    if 'colors_changed' not in st.session_state: 
        st.session_state.colors_changed = None

    if not st.session_state.files_uploaded:
        with st.form("file_upload_form"):
            counts_file = st.file_uploader("Upload the Counts file", type="csv")

            # Check if cell_annot_file has been uploaded
            if 'cell_annot_uploaded' not in st.session_state:
                st.session_state.cell_annot_uploaded = False
                st.session_state.df_cell_annot = None
            
            # first take file input
            if not st.session_state.cell_annot_uploaded:
                cell_annot_file = st.file_uploader("Upload the Cell Annotation file")
                if cell_annot_file is not None:
                    st.session_state.df_cell_annot = load_data(cell_annot_file)
                    st.session_state.cell_annot_uploaded = True
                    st.experimental_rerun()
            
            # Display dropdown menu with column names
            if st.session_state.cell_annot_uploaded:
                cell_annot_columns = st.session_state.df_cell_annot.columns.tolist()

                # Remove 'Unnamed: 0' if it exists and get the first column
                if 'Unnamed: 0' in cell_annot_columns:
                    cell_annot_columns.remove('Unnamed: 0')

                selected_column = st.selectbox("Select respective column for forming Sample Annotation table", 
                                               sorted(cell_annot_columns, key=str.lower),
                                               index=0)
                
                if selected_column:
                    # Extract unique categories from sample IDs
                    sample_ids = st.session_state.df_cell_annot[selected_column].astype(str).unique()
                    categories = set(sample_id.split('-')[0] for sample_id in sample_ids if '-' in sample_id)
                    
                    # Create selectboxes for case and control
                    case_category = st.selectbox("Select case category", options=sorted(list(categories)))
                    control_category = st.selectbox("Select control category", options=sorted(list(categories)))
            
            sample_annot_file = st.file_uploader("Upload the Sample Annotation file", type="txt")
            
            # If sample_annot_file is not uploaded, create it manually
            if not sample_annot_file and st.session_state.cell_annot_uploaded and selected_column:
                # Get unique sample IDs from the selected column
                unique_sample_ids = st.session_state.df_cell_annot[selected_column].astype(str).unique()
                
                df_sample_annot = pd.DataFrame({
                    'sample_ID': unique_sample_ids,
                    'case_or_control': [''] * len(unique_sample_ids),
                    'health_status': [''] * len(unique_sample_ids)
                })

                # Fill health_status and case_or_control columns
                for idx, sample_id in enumerate(df_sample_annot['sample_ID']):
                    if '-' in sample_id:
                        category = sample_id.split('-')[0]
                        df_sample_annot.loc[idx, 'health_status'] = category
                        if category == case_category:
                            df_sample_annot.loc[idx, 'case_or_control'] = 'case'
                        elif category == control_category:
                            df_sample_annot.loc[idx, 'case_or_control'] = 'control'
                                    
                st.session_state.df_sample_annot = df_sample_annot
            else:
                st.info("In case Sample Annotations to be generated automatically, after selecting the Cell Annotation file please follow the guidelines.")

            col1, col2, col3, col4 = st.columns([1, 1, 2, 1])
            with col1:
                submit_button = st.form_submit_button("Upload Files")
            with col2:
                default_files_button = st.form_submit_button("Demo Dataset")

            if submit_button:
                if st.session_state.cell_annot_uploaded and counts_file:
                    if sample_annot_file:
                        st.session_state.df_sample_annot = load_data(sample_annot_file)
                    # If sample_annot_file is not uploaded, use the manually created one
                    elif 'df_sample_annot' in st.session_state:
                        pass  # already df_sample_annot is manually created 
                    else:
                        st.error("Please create the Sample Annotation table manually or upload a file.")
                        return
                    
                    st.session_state.df_counts = load_data(counts_file)
                    st.session_state.df_database = load_data("./input_data/ready_database.txt")
                    st.session_state.files_uploaded = True
                    st.success("Each of the files uploaded successfully!")
                    time.sleep(3)
                    st.session_state.page = "Perform Analysis"
                    st.experimental_rerun()
                else:
                    st.error("Please upload all required files before submitting.")

            if default_files_button:    
                load_default_files()
                st.session_state.cell_annot_uploaded = True
                st.success("Demo dataset files loaded successfully!")
                time.sleep(3)
                st.session_state.page = "Perform Analysis"
                st.experimental_rerun()

    # once all the files are uploaded successfuly
    if st.session_state.files_uploaded:
        st.markdown("""
        <style>
            .stExpander {
                max-width: 100%;
                width: 100%;
            }
            .main .block-container {
                max-width: 80%;
                padding-top: 1rem;
                padding-right: 1rem;
                padding-left: 1rem;
                padding-bottom: 1rem;
            }
        </style>
        """, unsafe_allow_html=True)

        
        with st.expander("Data Files", expanded=True):
            database_tab, counts_tab, cell_annot_tab, sample_annot_tab = st.tabs(["Database", "Counts", "Cell Annotations", "Sample Annotations"])

            with database_tab:
                df_database = st.session_state.df_database

                # Multi-select for additional columns
                all_columns = df_database.columns.tolist()
                default_columns = ['Pair.Name', 'Ligand', 'Receptor']
                selected_columns = st.multiselect(
                    "Select additional columns to display:",
                    options=sorted([col for col in all_columns if col not in default_columns], key=str.lower),
                    default=[]
                )

                # Combine default and selected columns
                display_columns = default_columns + selected_columns

                # Filter DataFrame to include only selected columns
                df_display = df_database[display_columns]

                # Configure grid options
                gb = GridOptionsBuilder.from_dataframe(df_display)
                gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
                gb.configure_side_bar()
                gb.configure_selection('multiple', use_checkbox=True)
                gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                gb.configure_grid_options(columnSize="autoSize")
                gb.configure_grid_options(enableQuickFilter=True)
                gridOptions = gb.build()

                search = st.text_input("Search through table:", key="search_database")

                gridOptions['quickFilterText'] = search

                # Display the grid
                grid_response = AgGrid(
                    df_display,
                    gridOptions=gridOptions,
                    data_return_mode='AS_INPUT', 
                    update_mode=GridUpdateMode.MODEL_CHANGED,
                    fit_columns_on_grid_load=True,
                    theme='streamlit',
                    enable_enterprise_modules=True,
                    height=400,
                    width='100%',
                    reload_data=True
                )

                col1, col2, col3 = st.columns([1, 1, 1])
                with col1:
                    if st.button("Add New Entry"):
                        st.session_state.show_popup = "add"
                with col2:
                    if st.button("Delete Selected Rows"):
                        selected_rows = grid_response['selected_rows']
                        
                        if isinstance(selected_rows, pd.DataFrame) and not selected_rows.empty:
                            pair_names_to_delete = selected_rows['Pair.Name'].tolist()
                            
                            if pair_names_to_delete:
                                df_database = df_database[~df_database['Pair.Name'].isin(pair_names_to_delete)]
                                st.session_state.df_database = df_database
                                st.success(f"{len(pair_names_to_delete)} row(s) deleted successfully!")
                                time.sleep(3)
                                st.experimental_rerun()
                            else:
                                st.warning("No valid rows selected for deletion.")
                        else:
                            st.warning("No rows selected for deletion.")
                with col3:
                    if st.button("Delete Ligand/Receptor Family"):
                        st.session_state.show_popup = "delete_family"


                modal = Modal("Database Operations", key="database_modal")
                if st.session_state.get('show_popup') == "add":
                    with modal.container():
                        with st.form("new_entry_form"):
                            col1, col2 = st.columns(2)
                            with col1:
                                new_ligand = st.text_input("Ligand")
                            with col2:
                                new_receptor = st.text_input("Receptor")
                            
                            col3, col4, col5 = st.columns([1, 6, 1])
                            with col3:
                                submit_button = st.form_submit_button("Add")
                            with col5:
                                cancel_button = st.form_submit_button("Cancel")

                        if submit_button and new_ligand and new_receptor:
                            new_pair_name = f"{new_ligand}_{new_receptor}"
                            new_entry = pd.DataFrame({
                                'Pair.Name': [new_pair_name],
                                'Ligand': [new_ligand],
                                'Receptor': [new_receptor]
                            })

                            for col in df_database.columns:
                                if col not in new_entry.columns:
                                    new_entry[col] = ''

                            df_database = pd.concat([df_database, new_entry], ignore_index=True)
                            st.session_state.df_database = df_database
                            st.success(f"New entry added: {new_pair_name}")
                            time.sleep(3)
                            st.session_state.show_popup = False
                            modal.close()
                            st.experimental_rerun()

                        if cancel_button:
                            st.session_state.show_popup = False
                            modal.close()
                            st.experimental_rerun()

                elif st.session_state.get('show_popup') == "delete_family":
                    with modal.container():
                        with st.form("delete_family_form"):
                            delete_type = st.radio("Select type to delete", ["Ligand", "Receptor"])
                            delete_value = st.text_input("Enter the name:")
                            
                            col1, col2, col3 = st.columns([1, 6, 1])
                            with col1:
                                delete_button = st.form_submit_button("Delete")
                            with col3:
                                cancel_button = st.form_submit_button("Cancel")

                        if delete_button and delete_value:
                            if delete_type == "Ligand":
                                df_database = df_database[df_database['Ligand'] != delete_value]
                            else:
                                df_database = df_database[df_database['Receptor'] != delete_value]
                            
                            st.session_state.df_database = df_database
                            st.success(f"{delete_type} '{delete_value}' deleted from all instances.")
                            time.sleep(3)
                            st.session_state.show_popup = False
                            modal.close()
                            st.experimental_rerun()

                        if cancel_button:
                            st.session_state.show_popup = False
                            modal.close()
                            st.experimental_rerun()
            
            with counts_tab:
                df_counts = st.session_state.df_counts

                # Configure grid options
                gb = GridOptionsBuilder.from_dataframe(df_counts)
                gb.configure_pagination(paginationAutoPageSize=False)
                gb.configure_side_bar()
                gb.configure_selection('multiple')
                gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                gb.configure_grid_options(columnSize="autoSize")
                gb.configure_grid_options(enableQuickFilter=True)
                gridOptions = gb.build()

                search = st.text_input("Search through table:", key="search_counts")
                gridOptions['quickFilterText'] = search
                
                # Display the grid
                grid_response = AgGrid(
                    df_counts,
                    gridOptions=gridOptions,
                    data_return_mode='AS_INPUT', 
                    update_mode=GridUpdateMode.MODEL_CHANGED,
                    fit_columns_on_grid_load=True,
                    theme='streamlit',
                    enable_enterprise_modules=True,
                    height=400,
                    width='100%',
                    reload_data=True
                )

            with cell_annot_tab:
                df_cell_annot = st.session_state.df_cell_annot

                # Ensure 'sample_ID' is the first column
                if 'sample_ID' in df_cell_annot.columns:
                    cols = ['sample_ID'] + [col for col in df_cell_annot.columns if col != 'sample_ID']
                    df_cell_annot = df_cell_annot[cols]
                
                # Multi-select for additional columns
                all_columns = df_cell_annot.columns.tolist()
                default_columns = ['sample_ID', 'case_or_control', 'health_status']
                selected_columns = st.multiselect(
                    "Select additional columns to display:",
                    options=sorted([col for col in all_columns if col not in default_columns], key=str.lower),
                    default=[]
                )

                # Combine default and selected columns
                display_columns = default_columns + selected_columns

                # Filter DataFrame to include only selected columns
                df_display = df_cell_annot[display_columns]

                # Configure grid options
                gb = GridOptionsBuilder.from_dataframe(df_display)
                gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
                gb.configure_side_bar()
                gb.configure_selection('multiple', use_checkbox=True)
                gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                gb.configure_grid_options(columnSize="autoSize")
                gb.configure_grid_options(enableQuickFilter=True)
                gridOptions = gb.build()

                search = st.text_input("Search Across All Columns:", key="search_cell_annot")

                # Update gridOptions with the search text
                gridOptions['quickFilterText'] = search

                # Display the grid
                grid_response = AgGrid(
                    df_display,
                    gridOptions=gridOptions,
                    data_return_mode='AS_INPUT', 
                    update_mode=GridUpdateMode.MODEL_CHANGED,
                    fit_columns_on_grid_load=True,
                    theme='streamlit',
                    enable_enterprise_modules=True,
                    height=400,
                    width='100%',
                    reload_data=True
                )

            with sample_annot_tab:
                df_sample_annot = st.session_state.df_sample_annot

                # Ensure 'sample_ID' is the first column
                if 'sample_ID' in df_sample_annot.columns:
                    cols = ['sample_ID'] + [col for col in df_sample_annot.columns if col != 'sample_ID']
                    df_sample_annot = df_sample_annot[cols]
                
                # Multi-select for additional columns
                all_columns = df_sample_annot.columns.tolist()
                default_columns = ['sample_ID', 'case_or_control', 'health_status']
                selected_columns = st.multiselect(
                    "Select additional columns to display:",
                    options=sorted([col for col in all_columns if col not in default_columns], key=str.lower),
                    default=[]
                )

                # Combine default and selected columns
                display_columns = default_columns + selected_columns

                # Filter DataFrame to include only selected columns
                df_display = df_sample_annot[display_columns]

                # Configure grid options
                gb = GridOptionsBuilder.from_dataframe(df_display)
                gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
                gb.configure_side_bar()
                gb.configure_selection('multiple', use_checkbox=True)
                gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                gb.configure_grid_options(columnSize="autoSize")
                gb.configure_grid_options(enableQuickFilter=True)
                gridOptions = gb.build()

                search = st.text_input("Search Across All Columns:", key="search_sample_annot")

                # Update gridOptions with the search text
                gridOptions['quickFilterText'] = search

                # Display the grid
                grid_response = AgGrid(
                    df_display,
                    gridOptions=gridOptions,
                    data_return_mode='AS_INPUT', 
                    update_mode=GridUpdateMode.MODEL_CHANGED,
                    fit_columns_on_grid_load=True,
                    theme='streamlit',
                    enable_enterprise_modules=True,
                    height=400,
                    width='100%',
                    reload_data=True
                )
        
        with st.expander("Thresholds and Test Type", expanded=False):
            col1, col2, col3, col4, col5 = st.columns([3, 1, 3, 1, 3])
            with col1:
                threshold_celltype_size = st.slider("Cell type size", 0, 100, 6, help="`Cell type size` is a threshold for the minimum number of cells that a cell type should contain (in one sample). If the number of cells in the cell type of interest in a particular sample is less or equal to the `Cell type size`, then we consider this cell type as missing in this sample. This threshold affects the relative cell type abundance parameter (rho).")
                threshold_nr_active_cells = st.slider("Minimum number of active cells", 0, 100, 6, help="`Minimum number of active cells` is a threshold for the minimum number of active cells in a cell type (in the sample of interest). A cell is considered as active (for a specific gene), if it is expressing this gene above the `Threshold of expression`. If the number of active cells (for a specific gene) in a cell type is smaller or equal to the `Minimum number of active cells`, i.e. does not pass the threshold, then it is set to zero (in this sample). This threshold affects the relative active fraction (phi) parameter. The default value for the `Minimum number of active cells` is zero.")
                threshold_expr = st.slider("Threshold of expression", 0.0, 50.0, 0.1, step = 0.1, help="`Threshold of expression` is a threshold for an expression value of a gene in a cell. If an expression value af a gene in a cell is smaller or equal to the `Threshold of expression` value, it will be set to zero. This threshold affects the relative active fraction (phi) and the relative mean expression (p) parameters. The default value for the `Threshold of expression` is 0.05.")
                threshold_log10_cum_weight = st.slider("Log10 cumulative weight", 0.0, 50.0, 0.01, step = 0.01, help="For the quality check, we use three filters: the **interaction weight filter**, the **presence per cohort filter** and the **ligand/receptor expression filter**. An interaction is considered of good quality, if it passes all three filters. The **interaction weight filter** checks the log10 cumulative weight of the interaction. To pass this filter, the interaction needs to be greater than the `Log10 cumulative weight` threshold.")
                
            with col3:
                threshold_log10_meanexpr_per_condition = st.slider("Log10 mean expression per condition", 0.0, 50.0, 0.02, step = 0.01, help="The **ligand/receptor expression filter** checks the mean expression level of the ligand and the receptor of an interaction in the case and the control samples (separately).\n\n"
                                                                                                                                            "This filter uses a `Log10 mean expression per condition` threshold. \n"
                                                                                                                                            "For each interaction four values are checked:\n\n"
                                                                                                                                            "- log10 mean expression of the ligand in sending cells in control samples  \n"
                                                                                                                                            "- log10 mean expression of the receptor in receiving cells in control samples  \n"
                                                                                                                                            "- log10 mean expression of the ligand in sending cells in case samples  \n"
                                                                                                                                            "- log10 mean expression of the receptor in receiving cells in case samples.\n\n"
                                                                                                                                            "An interaction passes this filter if both its ligand and receptor pass the threshold either in control samples or in case samples or in both.")
                threshold_frac_samples_per_condition = st.slider("Fraction of samples per condition", 0.0, 50.0, 0.6, step = 0.1, help = "The **presence per cohort filter** checks the fraction of samples in which an interaction was detected (i.e. has a non-zero value) in the control cohort and in the case cohort. To pass this filter, an interaction needs to have a greater value than the `Fraction of samples per condition` threshold either in the control cohort or in the case cohort or in both.")
                threshold_log2FC = st.slider("Log2 fold change threshold", 1e-5, 1.0, 1.0, step=1e-5, help="For calculating statistically significant differential interactions between the cases and the controls, we need to define an adjusted p-value threshold and the log2 fold change threshold.")
                threshold_fdr = st.slider("False Discovery Rate (FDR) threshold", 0.0, 50.0, 0.1, step=0.1, help="Set the False Discovery Rate (FDR) threshold to control the expected proportion of false positives among significant results.")
            
            with col5:
                st.selectbox("Which test method?", ("Wilcoxon", "t-test"), index=1)


        if not st.session_state.analysis_performed:
            st.button("Perform Analysis", on_click=perform_analysis_callback)

        if st.session_state.analysis_performed:
            with st.expander("Calculate Communications", expanded=True):       
                with st.spinner("Running analysis..."):
                    if run_r_script("backend.R"):
                        st.success("Analysis completed successfully!")
                        
                        try:
                            interactions_df = pd.read_csv("interactions_QC.txt", sep="\t")
                            
                            st.subheader("Interaction Results")
                            
                            # Configure grid options for interactions table
                            gb = GridOptionsBuilder.from_dataframe(interactions_df)
                            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
                            gb.configure_side_bar()
                            gb.configure_selection('single')
                            gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                            gb.configure_grid_options(domLayout='normal')
                            gridOptions = gb.build()

                            # Display the interactions grid
                            AgGrid(
                                interactions_df,
                                gridOptions=gridOptions,
                                data_return_mode='AS_INPUT', 
                                update_mode=GridUpdateMode.MODEL_CHANGED,
                                fit_columns_on_grid_load=False,
                                theme='streamlit',
                                enable_enterprise_modules=True,
                                height=400,
                                width='100%',
                                reload_data=True,
                                columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS
                            )

                            st.subheader("Filter Interactions")

                            comm_result = prepare_comm_result(interactions_df)

                            filtered_result, fig_cumW, fig_fracSamp, fig_meanLigRec = filter_interactions(comm_result)

                            #st.plotly_chart(fig_cumW)
                            #st.plotly_chart(fig_fracSamp)
                            #st.plotly_chart(fig_meanLigRec)

                            image_dir = "./plots/"
                            col1, col2 = st.columns(2)
                            with col1:
                                load_and_display_image(os.path.join(image_dir, "above_plot1.png"))
                                load_and_display_image(os.path.join(image_dir, "meanlig-vs-meanrec.png"))
                            with col2:
                                load_and_display_image(os.path.join(image_dir, "above_plot2.png"))
                            
                            # Display filtered results
                            st.subheader("Filtered Interactions")
                            filtered_interactions = filtered_result['anno_interactions'][filtered_result['anno_interactions']['passed_QC_filter']]
                            
                            gb = GridOptionsBuilder.from_dataframe(filtered_interactions)
                            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
                            gb.configure_side_bar()
                            gb.configure_selection('single')
                            gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                            gb.configure_grid_options(domLayout='normal')
                            gridOptions = gb.build()

                            AgGrid(
                                filtered_interactions,
                                gridOptions=gridOptions,
                                data_return_mode='AS_INPUT', 
                                update_mode=GridUpdateMode.MODEL_CHANGED,
                                fit_columns_on_grid_load=False,
                                theme='streamlit',
                                enable_enterprise_modules=True,
                                height=400,
                                width='100%',
                                reload_data=True,
                                columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS
                            )

                            diff_exp_interactions_df = pd.read_csv("diff_exp_interactions.txt", sep="\t")

                            st.subheader("Differential Interactions")
                            
                            # Configure grid options for interactions table
                            gb = GridOptionsBuilder.from_dataframe(diff_exp_interactions_df)
                            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
                            gb.configure_side_bar()
                            gb.configure_selection('single')
                            gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc="sum", editable=True)
                            gb.configure_grid_options(domLayout='normal')
                            gridOptions = gb.build()

                            # Display the interactions grid
                            AgGrid(
                                diff_exp_interactions_df,
                                gridOptions=gridOptions,
                                data_return_mode='AS_INPUT', 
                                update_mode=GridUpdateMode.MODEL_CHANGED,
                                fit_columns_on_grid_load=False,
                                theme='streamlit',
                                enable_enterprise_modules=True,
                                height=400,
                                width='100%',
                                reload_data=True,
                                columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS
                            )

                        except FileNotFoundError:
                            st.error("interactions.txt file not found. Make sure the R script generates this file.")
                        except pd.errors.EmptyDataError:
                            st.warning("The interactions.txt file is empty.")
                        except Exception as e:
                            st.error(f"An error occurred while reading the interactions file: {e}")

            with st.expander("Visualizations", expanded=True):
                
                #aml_data = pd.read_csv("visualizations/AML_data.txt", sep="\t")
                #healthy_data = pd.read_csv("visualizations/Healthy_data.txt", sep="\t")
                
                #fig_aml, fig_healthy = plot_aml_healthy_visual(aml_data, healthy_data)
                #st.plotly_chart(fig_aml)
                #st.plotly_chart(fig_healthy)

                if st.session_state.colors_changed == None:
                    run_r_script("visualization.R")
                display_visualizations()

                
def main():
    #st.set_page_config(initial_sidebar_state="collapsed")
    
    if 'page' not in st.session_state:
        st.session_state.page = "Home"

    st.sidebar.title("Navigation to")
    pages = ["Home", "Perform Analysis", "License"]
    page = st.sidebar.radio("", pages, index=pages.index(st.session_state.page))

    if page != st.session_state.page:
        st.session_state.page = page
        st.rerun()

    if st.session_state.page == "Home":
        show_home_page()
    elif st.session_state.page == "Perform Analysis":
        perform_analysis()
    elif st.session_state.page == "License":
        show_license_page()


if __name__ == "__main__":
    main()