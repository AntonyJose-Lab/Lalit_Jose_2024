import pandas as pd
import numpy as np

import pickle
import networkx as nx
import plotly.graph_objects as go

import dash
from dash import dcc, html
from dash.dependencies import Input, Output

import community

# Distance matrix

dist_matrix = pd.read_csv("dist_matrix.csv", index_col=0)


# DASH APP

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div([
    html.H1("Gene Network Visualization"),
    
    html.Div([
        html.Label("Enter list of genes to be highlighted on the network: "),
        html.Br(),
        html.Br(),
        dcc.Textarea(
            id='gene-input',
            placeholder='Separate genes with commas...',
            style={'width': '50%'}
        ),
        html.Br(),
        html.Br(),
        html.Label("Enter threshold as a real number between 0 and 1: "),
        html.Br(),
        html.Br(),
        dcc.Input(
            id='threshold-input',
            value='0.75', 
            type='number',
            step='any', 
        ),
        dcc.Store(id='stored-graph-data'),
        html.Br(),
        html.Br(),
        html.Button('Submit', id='submit-button', n_clicks=0),
        html.Br(),
        html.Br()
    ]),
    dcc.Graph(id='network-graph', figure={
    'data': [],
    'layout': go.Layout(
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, 
                   visible=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, 
                   visible=False)
    )
}),
    html.Br(),
    html.Div(id='feedback'),
    html.Br()
])

def generate_contrasting_colors(n):
    colors = []
    while len(colors) < n:
        r, g, b = np.random.rand(3)
        if not (g > r and g > b) and f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' not in colors:
            colors.append(f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}')
            
    return colors


def create_and_save_community_df(partition):

    # Reverse the partition dictionary to group nodes by community
    reversed_partition = {}
    for node, comm_id in partition.items():
        if comm_id not in reversed_partition:
            reversed_partition[comm_id] = []
        reversed_partition[comm_id].append(node)

    # Filter out the single-node communities
    filtered_partition = {k: v for k, v in reversed_partition.items() if len(v) > 1}

    # Create the DataFrame
    community_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in filtered_partition.items()]))

    # Save the DataFrame to CSV
    community_df.to_csv("communities.csv", index=False)

def save_graph_positions_and_clusters(G, pos, partition, file_path="graph_data.pickle"):
    with open(file_path, 'wb') as f:
        pickle.dump((G, pos, partition), f)


def consistent_girvan_newman_partition(G, num_iterations=10):
    all_communities = []  # To store sets of nodes for each community in each iteration

    for _ in range(num_iterations):
        communities_generator = nx.community.girvan_newman(G)
        top_level_communities = next(communities_generator)

        community_groups = {}
        for idx, community in enumerate(top_level_communities):
            for node in community:
                if idx not in community_groups:
                    community_groups[idx] = set()
                community_groups[idx].add(node)
        
        # Store the sets of nodes for this iteration
        all_communities.append(list(community_groups.values()))

    # Determine nodes that consistently appear together
    consistent_nodes = set(G.nodes())
    for node in G.nodes():
        communities_node_belongs_to = [comm for community_list in all_communities for comm in community_list if node in comm]
        
        if not all(communities_node_belongs_to[0] == community for community in communities_node_belongs_to):
            consistent_nodes.remove(node)

    # Generate the final partition using Girvan-Newman
    communities_generator = nx.community.girvan_newman(G)
    top_level_communities = next(communities_generator)
    final_partition = {}
    for idx, community in enumerate(top_level_communities):
        for node in community:
            if node in consistent_nodes:
                final_partition[node] = idx
            else:
                final_partition[node] = -1  # Inconsistent nodes

    return final_partition




@app.callback(
    [Output('network-graph', 'figure'),
     Output('feedback', 'children'),
     Output('stored-graph-data', 'data')],
    [Input('submit-button', 'n_clicks')],
    [dash.dependencies.State('gene-input', 'value'),
     dash.dependencies.State('threshold-input', 'value'),
     dash.dependencies.State('stored-graph-data', 'data')]
)

def update_graph(n_clicks, gene_input, threshold_value, stored_graph_data):
    
    threshold = float(threshold_value)

    # Initialize stored_graph_data if None or empty
    if stored_graph_data is None or not stored_graph_data:
        stored_graph_data = {}

    # Check if the graph needs to be regenerated
    regenerate_graph = stored_graph_data.get('prev_threshold') != threshold

    if regenerate_graph:  # Only regenerate graph and recolor if the threshold changes.
        adjacency_matrix = (dist_matrix < threshold).astype(int)
        np.fill_diagonal(adjacency_matrix.values, 0)
        G = nx.from_pandas_adjacency(adjacency_matrix)
        pos = nx.spring_layout(G, k=0.15)

        # Community detection and coloring
        partition = consistent_girvan_newman_partition(G)
        create_and_save_community_df(partition)
        communities = set(partition.values())
        colors = generate_contrasting_colors(len(communities))
        community_colors = {comm_id: colors[i] for i, comm_id in enumerate(communities)}
        

        # Store the generated graph data, position data, community partition, and community colors for later use
        stored_graph_data['adj_list'] = nx.to_dict_of_lists(G)
        stored_graph_data['pos'] = pos
        stored_graph_data['partition'] = partition
        stored_graph_data['community_colors'] = community_colors
        stored_graph_data['prev_threshold'] = threshold

    else:  # If the threshold hasn't changed, retrieve the stored graph data.
        G = nx.from_dict_of_lists(stored_graph_data['adj_list'])
        pos = stored_graph_data['pos']
        partition = stored_graph_data['partition']
        community_colors = stored_graph_data['community_colors']

    # Highlight the nodes based on the input list of genes. No regeneration of graph or recoloring.
    highlighted_genes = [gene.strip().lower() for gene in gene_input.split(",")] if gene_input else []
    not_in_graph = [gene for gene in highlighted_genes if gene.lower() not in [g.lower() for g in G.nodes()]]
    feedback_msg = ', '.join(not_in_graph) + " are not in the graph." if not_in_graph else ""

    for node in G.nodes():
        if G.degree(node) == 0:  # Check for unconnected nodes
            community_colors[partition[node]] = '#000000'  # Color them black
    
    # Assigning black color to inconsistent nodes (those marked with -1)
    community_colors[-1] = '#808080' 
    community_colors[-2] = '#000000'

    for node, comm_id in partition.items():
        if comm_id not in community_colors:
            # Generate non-green colors
            r, b = (int(255 * np.random.rand()), int(255 * np.random.rand()))
            g = int(min(r, b) * 0.3)
            community_colors[comm_id] = f'#{r:02x}{g:02x}{b:02x}'
            
    node_colors_regular = []        
    node_colors_highlighted = []

    node_x_regular = []
    node_y_regular = []
    node_labels_regular = []

    node_x_highlighted = []
    node_y_highlighted = []
    node_labels_highlighted = []

    for k in G.nodes():
        neighbors = list(G.neighbors(k))
        neighbors_sorted_by_distance = sorted(neighbors, key=lambda x: dist_matrix.at[k, x])
        top_neighbors = neighbors_sorted_by_distance[:10]
        hover_info = f"Gene: {k}<br>Cluster ID: {partition[k]}<br><br>Connected to (Top 10):<br>{'<br>'.join(top_neighbors)}"
        
        # If the current node is one of the highlighted genes, it is colored green
        if k.lower() in highlighted_genes:
            node_x_highlighted.append(pos[k][0])
            node_y_highlighted.append(pos[k][1])
            node_labels_highlighted.append(hover_info)
            node_colors_highlighted.append("green")
        else:
            # If not highlighted, we determine and assign the community color for the node.
            color = community_colors[partition[k]]
            node_x_regular.append(pos[k][0])
            node_y_regular.append(pos[k][1])
            node_labels_regular.append(hover_info)
            node_colors_regular.append(color)

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')
    
    node_trace_regular = go.Scatter(
        x=node_x_regular, y=node_y_regular, mode='markers', hoverinfo='text', marker=dict(showscale=False, size=10, line_width=2, color=node_colors_regular), text=node_labels_regular)
    
    node_trace_highlighted = go.Scatter(
        x=node_x_highlighted, y=node_y_highlighted, mode='markers', hoverinfo='text', marker=dict(size=10, line_width=2, color=node_colors_highlighted), text=node_labels_highlighted)
    
    return go.Figure(data=[edge_trace, node_trace_regular, node_trace_highlighted],
                    layout=go.Layout(title='Historical mutual information of genes in RNA silencing',
                                    titlefont_size=16,
                                    showlegend=False,
                                    hovermode='closest',
                                    margin=dict(b=0, l=0, r=0, t=40),
                                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))), \
        feedback_msg, stored_graph_data

if __name__ == '__main__':
    app.run_server(debug=False)
