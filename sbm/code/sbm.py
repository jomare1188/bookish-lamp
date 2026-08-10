import logging
import numpy as np
from collections import defaultdict
import argparse
import sys
from random import choice
from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

import pandas as pd

# Logging configuration
FORMAT = "%(message)s"
# format="[%(levelname)s]: "
logging.basicConfig(
    format=FORMAT,
    level="INFO",
    handlers=[RichHandler(show_time=False, show_path=False, markup=True)],
)

logtext = logging.getLogger("rich")
logtext.setLevel(20)

console = Console()

try:
    import graph_tool.all as gt
    from graph_tool.all import *
except ImportError:
    logtext.error(
        "graph-tool is not installed. Please install it to use SBM functionalities.\nCheck if the conda environment is activated.\nSee https://graph-tool.skewed.de/installation.html"
    )
    sys.exit(1)

def _make_progress() -> Progress:
    """
    Create a rich Progress object with a random spinner and various progress columns.
    Returns:
        Progress: A rich Progress object for displaying progress bars in the console.
    """
    spinners = [
        "aesthetic",
        "shark",
        "dots",
        "line",
        "bouncingBall",
        "moon",
        "earth",
        "monkey",
        "runner",
        "pong",
        "weather",
        "clock",
    ]
    return Progress(
        SpinnerColumn(choice(spinners)),
        TaskProgressColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        transient=True,
    )


class BioLayerNetLogger:
    """
    A simple logger class for logging messages with different severity levels (info, warning, error) using the rich library.
    """

    @staticmethod
    def info(message: str) -> None:
        """Log an info message."""
        logtext.info(message)

    @staticmethod
    def warning(message: str) -> None:
        """Log a warning message."""
        logtext.warning(message)

    @staticmethod
    def error(message: str) -> None:
        """Log an error message."""
        logtext.error(message)


logger = BioLayerNetLogger()


def check_weighted_model(
    graph: gt.Graph, weight_prop: gt.PropertyMap
) -> tuple[gt.NestedBlockState, gt.NestedBlockState]:
    """
    Check the model fit with edge weights (covariates) by comparing the description lengths of a standard SBM and a weighted SBM.
    In the weighted SBM, the edge weights (covariates) are taken into account, while in the standard SBM, they are ignored.
    The model with the lower description length is preferred.

    Args:
        graph (gt.Graph): The input graph-tool graph.
        weight_prop (gt.PropertyMap): Edge property map indicating the weights (covariates) of each edge.
    Returns:
        tuple[gt.NestedBlockState, gt.NestedBlockState]: A tuple containing the preferred nested block state (either standard or weighted) and the standard state for comparison.
    """
    logger.info(
        "Checking model fit with edge weights (covariates)... This usually takes longer than the standard model. Take a while... Get a coffee :)"
    )

    # PTBR
    # Modelo Padrão (Topológico Puro)
    # Serve como base para ver como a rede se agrupa apenas por conexões (quem liga com quem)
    # EN
    # Standard Model (Pure Topology)
    # Serves as a baseline to see how the network groups based solely on connections (who connects to whom)
    logger.info("Running standard (unweighted) model...")
    with _make_progress() as p:
        p.add_task(description="[cyan]Inferring standard SBM...", total=None)
        state_std = gt.minimize_nested_blockmodel_dl(
            graph, state_args=dict(deg_corr=True)
        )
    dl_std = state_std.entropy()
    logger.info(f"Standard model Description Length (DL): {dl_std:.2f}")

    # PTBR
    # Modelo Ponderado (Weighted SBM)
    # Utiliza 'recs' (covariáveis) em vez de 'layers'.
    # Assumimos distribuição Poisson pois seus dados são contagens (0, 1, 2...)
    # EN
    # Weighted Model (Weighted SBM)
    # Uses 'recs' (covariates) instead of 'layers'.
    # We assume a Poisson distribution since your data are counts (0, 1, 2...)
    logger.info("Running weighted model...")
    with _make_progress() as p:
        p.add_task(description="[cyan]Inferring weighted SBM...", total=None)
        state_weighted = gt.minimize_nested_blockmodel_dl(
            graph,
            #base_state=gt.WeightedBlockState,
            state_args=dict(
                deg_corr=True,
                recs=[
                    weight_prop
                ],  # Propriedade dos pesos. Correlecao de Pearson com valor absoluto > 0.8 e depois MinMax (= epsilon para evitar 0)
                rec_types=[
                    "real-exponential"
                ],  # Distribuição adequada para os dados. Verificar histograma dos pesos.
            ),
        )

    dl_weighted = state_weighted.entropy()
    logger.info(f"Weighted model Description Length (DL): {dl_weighted:.2f}")

    # PTBR
    # NOTA CRÍTICA SOBRE COMPARAÇÃO DE DL:
    # Ao contrário do modelo Layered vs Standard onde a estrutura de dados era similar,
    # o modelo Ponderado tem MUITO mais dados para descrever (os valores das arestas).
    # Por isso, o DL_weighted quase sempre será MAIOR que o DL_std.
    # A comparação direta 'if dl_weighted < dl_std' não é justa matematicamente aqui.
    # O teste real seria: "O modelo ponderado encontrou uma partição diferente?"
    # Caso você queira usar esses atributos, o estado ponderado é o que contém a informação rica.
    # EN
    # CRITICAL NOTE ON DL COMPARISON:
    # Unlike the Layered vs Standard model where the data structure was similar,
    # the Weighted model has MUCH more data to describe (the edge values).
    # Therefore, the DL_weighted will almost always be GREATER than the DL_std.
    # The direct comparison 'if dl_weighted < dl_std' is not mathematically fair here.
    # The real test would be: "Did the weighted model find a different partition?"
    # If you want to use these attributes, the weighted state is the one that contains the rich information.

    logger.info(
        "[green]Weighted model inference complete. Using this state implies that blocks are grouped "
        "by both connectivity pattern AND interaction strength (BioGRID counts).[/]"
    )
    # PTBR
    # Retornamos o ponderado pois ele contém a informação biológica do grafo
    # EN
    # We return the weighted state because it contains the biological information of the graph
    return state_weighted, state_std


def sbm_annealing(
    state: gt.NestedBlockState,
    anneal_niter: int = 100,
    equilibrate_wait: int = 100,
    marginal_niter: int = 100,
    mcmc_niter: int = 100,
    verbose: bool = False,
) -> gt.NestedBlockState:
    """
    Executes the SBM annealing process, including MCMC annealing, equilibration, and marginal sampling.
    This process is computationally intensive and may take a while, so it's recommended to have a coffee ready!
    This function executes the following steps:
    1. MCMC Annealing: Gradually cools the system to find a good partitioning of the graph.
    2. MCMC Equilibration: Further refines the partitioning by allowing the system to reach equilibrium.
    3. MCMC Marginal Sampling: Collects samples of partitions to compute a consensus partition.

    Args:
        state (gt.NestedBlockState): The initial nested block state.
        anneal_niter (int): Number of iterations for the annealing process.
        equilibrate_wait (int): Number of iterations to wait for equilibration.
        marginal_niter (int): Number of iterations for marginal sampling.
        mcmc_niter (int): Number of MCMC iterations per step.
        verbose (bool): If True, enables verbose output during the process.
    Returns:
        gt.NestedBlockState: The final nested block state after annealing, equilibration, and marginal sampling.
    """
    logger.info(
        "☕ [cyan]Starting SBM annealing process… This may take a while. Perfect time to refill your coffee![/]"
    )

    state = state.copy()
    s1_entropy_initial = state.entropy()
    logger.info(f"Initial state entropy: {s1_entropy_initial:.2f}")
    logger.info("🔥 Beginning MCMC annealing…")
    with _make_progress() as p:
        p.add_task(description="[cyan]Running MCMC annealing...", total=None)
        gt.mcmc_anneal(
            state,
            beta_range=(1.0, 10.0),
            niter=anneal_niter,
            mcmc_equilibrate_args=dict(force_niter=mcmc_niter),
            verbose=verbose,
        )
    s2_entropy = state.entropy()
    logger.info(f"Post-annealing state entropy: {s2_entropy:.2f}")
    logger.info(
        f"Entropy improvement after annealing: {s1_entropy_initial - s2_entropy:.2f}"
    )

    logger.info("🌀 Starting MCMC equilibration… Time for another sip of coffee ☕")

    s1_entropy = state.entropy()
    with _make_progress() as p:
        p.add_task(description="[cyan]Running MCMC equilibration...", total=None)
        gt.mcmc_equilibrate(
            state,
            wait=equilibrate_wait,
            mcmc_args=dict(niter=mcmc_niter),
            verbose=verbose,
        )
    s2_entropy = state.entropy()
    logger.info(f"Post-equilibration state entropy: {s2_entropy:.2f}")
    logger.info(
        f"[green]Entropy improvement after equilibration: {s1_entropy - s2_entropy:.2f}"
    )

    logger.info(
        "🎲 Starting MCMC marginal sampling… This one's long — maybe grab a snack 🍪"
    )

    partition_samples = []

    def _collect_partitions(s):
        partition_samples.append(s.get_bs())

    with _make_progress() as p:
        p.add_task(description="[cyan]Running MCMC marginal sampling...", total=None)
        gt.mcmc_equilibrate(
            state,
            force_niter=marginal_niter,
            mcmc_args=dict(niter=mcmc_niter),
            callback=_collect_partitions,
            verbose=verbose,
        )

    logger.info(
        f"Collected {len(partition_samples)} partition samples during marginal sampling."
    )

    logger.info("🧠 Computing consensus partition… Hang tight!")
    pmode = gt.PartitionModeState(partition_samples, nested=True, converge=True)
    bs_consensus = pmode.get_max_nested()

    state = state.copy(bs=bs_consensus)
    s_final_entropy = state.entropy()
    logger.info(f"Final consensus state entropy: {s_final_entropy:.2f}")
    logger.info(
        f"Total entropy improvement: {s1_entropy_initial - s_final_entropy:.2f}"
    )
    return state


def get_uncertainty_reconstruction(
    graph: gt.Graph,
    initial_state: gt.NestedBlockState,
    equilibrate_wait: int = 100,
    mcmc_niter: int = 100,
    verbose: bool = True,
) -> gt.Graph:
    """
    Performs network reconstruction to quantify edge uncertainty and predict missing links.

    Args:
        graph: The input graph (BioGRID/KEGG).
        initial_state: (Optional) The high-quality partition from 'sbm_annealing'.
                             Passing this warms-up the reconstruction significantly.
        equilibrate_wait: Number of iterations to wait for equilibration.
        mcmc_niter: Number of MCMC iterations per step.
        verbose: If True, enables verbose output during the process.
    Returns:
        gt.Graph: A new graph containing the marginal probabilities of edges existing,
                  with the property 'eprob' representing the posterior probability of each edge.
    """

    logger.info(
        "🕵️ [cyan]Initializing Network Reconstruction (MeasuredBlockState)...[/]"
    )

    # We assume a "Simple Noise Model":
    # The graph structure itself is the observation 'x'.
    # We assume every pair was measured once (n=1), which is the default for MeasuredBlockState
    # if 'n' is not provided.

    # We define the "Simple Noise Model" for existing edges:
    # For every edge present in the graph, we claim we measured it ONCE (n=1)
    # and it was POSITIVE (x=1).
    n_prop = graph.new_edge_property("long")  # Must be 'long' (int64_t)
    x_prop = graph.new_edge_property("long")
    n_prop.a = 1
    x_prop.a = 1

    # Instead of passing the full 'initial_state' object (which risks graph mismatches),
    # we extract the partition hierarchy (b) as a list of arrays.
    init_b = None
    if initial_state:
        # Extract the block labels from each level of the hierarchy
        init_b = [lev.b.a.copy() for lev in initial_state.get_levels()]
        logger.info(
            "    Using optimized partition from annealing to warm up inference."
        )

    state = gt.MeasuredBlockState(
        graph,
        n=n_prop,
        x=x_prop,
        nested=True,
        state_args=dict(deg_corr=True, bs=init_b),
    )

    # 1. Equilibrate the reconstruction state
    # This infers the false-positive (q) and false-negative (p) rates
    # compatible with the SBM structure.
    logger.info("⚖️  Equilibrating reconstruction MCMC (inferring error rates p, q)...")
    with _make_progress() as p:
        p.add_task(
            description="[cyan]Equilibrating reconstruction state...", total=None
        )
        gt.mcmc_equilibrate(
            state,
            wait=equilibrate_wait,
            mcmc_args=dict(niter=mcmc_niter),
            verbose=verbose,
        )

    # 2. Collect Marginals
    # We sample the posterior to find the probability of every edge existing.
    logger.info(
        f"📊 Collecting marginal probabilities over {equilibrate_wait} iterations..."
    )

    u = None  # This will hold our marginal graph

    def _collector(s):
        """
        Callback function to collect marginals.
        """
        nonlocal u
        u = s.collect_marginal(u)

    with _make_progress() as p:
        p.add_task(description="[cyan]Sampling marginals...", total=None)
        gt.mcmc_equilibrate(
            state,
            force_niter=equilibrate_wait,
            mcmc_args=dict(niter=mcmc_niter),
            callback=_collector,
            verbose=verbose,
        )

    # Finalize the marginal graph (converts counts to probabilities)
    # Note: 'u' is already the marginal graph object updated in-place,
    # but we ensure we return the latest state.
    logger.info("✅ Reconstruction complete. Marginal graph generated.")
    return u


def integrate_reconstruction_results(
    graph: gt.Graph,
    marginal_graph: gt.Graph,
    block_state: gt.NestedBlockState = None,
    prob_prop: str = "eprob",
    prediction_threshold: float = 0.01,  # Lower threshold to see more potential links in DF
) -> tuple[gt.Graph, pd.DataFrame]:
    """
    Integrates reconstruction results into the graph and generates a comprehensive report.

    Outputs:
    1. Modifies 'graph' in-place: Adds 'posterior_prob' and 'reliability' properties.
    2. Returns 'df_results': A unified DataFrame of EXISTING and PREDICTED edges
       with block context (Source_Block, Target_Block, etc.).

    Args:
        graph: The original biological network.
        marginal_graph: Output from run_uncertainty_reconstruction.
        block_state: The annealed state (for block IDs).
        prob_prop: The edge property name in 'marginal_graph' that contains the posterior probabilities.
        prediction_threshold: Min probability to include a NEW PREDICTION in the DataFrame.
                              (Existing edges are ALWAYS included, even if prob=0).
    Returns:
        graph: The original graph with added edge properties.
        df_results: A pandas DataFrame summarizing existing and predicted edges with their probabilities and block context.
    """
    logger.info(" Processing reconstruction results (Unified Report)...")

    # Setup Graph Properties
    ep_prob = graph.new_edge_property("double")
    ep_reliability = graph.new_edge_property("string")
    graph.ep["posterior_prob"] = ep_prob
    graph.ep["reliability"] = ep_reliability

    # Access marginal probabilities
    if prob_prop not in marginal_graph.ep:
        raise ValueError(f"Property '{prob_prop}' not found in marginal graph.")
    mgr_prob = marginal_graph.ep[prob_prop]

    # Node names map
    v_name = graph.vp["name"] if "name" in graph.vp else graph.vertex_index

    # Block ID Helper
    def get_block_info(node_idx):
        if block_state:
            # Get Level 0 (lowest level) block ID
            return block_state.get_levels()[0].get_blocks()[node_idx]
        return -1

    # Data container for DataFrame
    all_edges_data = []

    # Helper to create a row for the DataFrame
    def create_row(src, tgt, prob, edge_type, status) -> dict:
        """
        Creates a dictionary representing a row in the results DataFrame for a given edge.

        Args:
            src (int): Source vertex index.
            tgt (int): Target vertex index.
            prob (float): Posterior probability of the edge.
            edge_type (str): "Existing" or "Predicted".
            status (str): Reliability status (e.g., "High Confidence", "Spurious").
        Returns:
            dict: A dictionary representing the edge information for the DataFrame.
        """
        row = {
            "Source_Name": v_name[src],
            "Target_Name": v_name[tgt],
            "Type": edge_type,  # "Existing" or "Predicted"
            "Status": status,  # "High Confidence", "Spurious", etc.
            "Posterior_Probability": prob,
        }

        # Add Block Context
        if block_state:
            b_src = get_block_info(src)
            b_tgt = get_block_info(tgt)
            row["Source_Block"] = b_src
            row["Target_Block"] = b_tgt
            # Useful for pivot tables: "10 <-> 15"
            row["Block_Interaction"] = f"{min(b_src, b_tgt)} <-> {max(b_src, b_tgt)}"

        return row

    # Track processed edges to avoid duplicates
    # Format: (min(u,v), max(u,v)) assuming undirected behavior for tracking
    # If directed, just (u,v). Assuming Directed for biological networks usually.
    processed_edges = set()

    # Iterate over MARGINAL Graph (High Prob Existing + Predictions)
    for e_marg in marginal_graph.edges():
        src, tgt = int(e_marg.source()), int(e_marg.target())
        prob = mgr_prob[e_marg]

        # Mark as processed
        processed_edges.add((src, tgt))

        # Check if edge exists in ORIGINAL
        e_orig = graph.edge(src, tgt)

        if e_orig is not None:
            # EXISTING EDGE
            ep_prob[e_orig] = prob

            if prob >= 0.9:
                status = "High Confidence"
            elif prob >= 0.5:
                status = "Uncertain"
            else:
                status = "Likely Spurious"

            ep_reliability[e_orig] = status

            # Always add existing edges to DataFrame
            all_edges_data.append(create_row(src, tgt, prob, "Existing", status))

        else:
            # NEW PREDICTION
            if prob >= prediction_threshold:
                all_edges_data.append(
                    create_row(src, tgt, prob, "Predicted", "New Interaction")
                )

    # Iterate over ORIGINAL Graph (To catch Zero-Prob/Spurious edges)
    # Some existing edges might be so unlikely they never appeared in the MCMC samples (prob=0)
    for e in graph.edges():
        src, tgt = int(e.source()), int(e.target())

        if (src, tgt) not in processed_edges:
            # This edge exists in input but was NOT in marginals -> Prob = 0.0
            prob = 0.0
            status = "Spurious (Prob ~ 0)"

            # Update Graph
            ep_prob[e] = prob
            ep_reliability[e] = status

            # Add to DataFrame
            all_edges_data.append(create_row(src, tgt, prob, "Existing", status))

    # Finalize DataFrame
    df_results = pd.DataFrame(all_edges_data)

    if not df_results.empty:
        # Reorder columns for readability
        cols = ["Source_Name", "Target_Name", "Type", "Status", "Posterior_Probability"]
        if block_state:
            cols += ["Source_Block", "Target_Block", "Block_Interaction"]

        # Sort: High probability first
        df_results = df_results[cols].sort_values(
            by="Posterior_Probability", ascending=False
        )

    logger.info("✅ Integrated results.")
    logger.info(
        f"DataFrame contains {len(df_results)} rows (Existing + Predicted > {prediction_threshold})"
    )

    return graph, df_results


def get_block_attrib(
    graph: gt.Graph, unique_blocks: list, block_array: gt.PropertyArray
) -> defaultdict:
    """
    Get the block statistics for each block in the graph, including size, internal edges, external edges, degrees, density, and q_r.

    Args:
        graph (gt.Graph): The input graph-tool graph.
        unique_blocks (list): List of unique block IDs.
        block_array (gt.PropertyArray): Array mapping each vertex to its block ID.
    Returns:
        defaultdict: A dictionary containing block statistics for each block ID.
    """
    block_stats = defaultdict(dict)
    E = graph.num_edges()
    directed = graph.is_directed()
    ba = np.asarray(block_array)

    # Pre-compute degree arrays once (C++ backed, O(V))
    if directed:
        in_deg = np.asarray(graph.get_in_degrees(graph.get_vertices()))
        out_deg = np.asarray(graph.get_out_degrees(graph.get_vertices()))
    else:
        out_deg = np.asarray(graph.get_out_degrees(graph.get_vertices()))
        in_deg = out_deg  # Symmetric for undirected

    # Get full edge list once as NumPy array (C++ backed, O(E))
    edges = graph.get_edges()  # shape (E, 2+)
    src_blocks = ba[edges[:, 0].astype(int)]
    tgt_blocks = ba[edges[:, 1].astype(int)]

    for block_id in unique_blocks:
        mask = ba == block_id
        nodes_in_block = np.where(mask)[0]
        n = len(nodes_in_block)

        # Internal edges: both endpoints in this block (vectorised O(E))
        internal_mask = (src_blocks == block_id) & (tgt_blocks == block_id)
        e_rr = int(np.sum(internal_mask))
        if not directed:
            e_rr = e_rr // 2  # Each undirected edge counted twice

        # Block degree sums (vectorised O(V))
        k_in = int(np.sum(in_deg[mask]))
        k_out = int(np.sum(out_deg[mask]))

        # Expected edges under configuration-model null
        if E > 0:
            if directed:
                expected_edges = (k_in * k_out) / E
            else:
                # undirected null: (sum of block degrees)^2 / (2 * 2E)  -> k_r^2/(2E)
                k_r = k_in
                expected_edges = (k_r * k_r) / (2 * E)
            # q_r for directed graph
            # based in the Fig4 https://arxiv.org/pdf/2006.14493
            # TODO: Check the generalization for undirected graphs
            q_r = (e_rr - expected_edges) / E
        else:
            q_r = 0

        if not (-1.0001 <= q_r <= 1.0001):
            logger.warning(
                f"q_r out of [-1, 1] for block {block_id}: {q_r:.4f}. "
                "This usually means the graph directedness was mis-detected.\n"
                " Please, open an issue in code repository."
            )

        if directed:
            external_edges = k_out - e_rr
            max_possible_edges = n * (n - 1)
        else:
            external_edges = k_out - 2 * e_rr
            max_possible_edges = n * (n - 1) / 2

        density = e_rr / max_possible_edges if max_possible_edges > 0 else 0

        block_stats[block_id] = {
            "block_id": int(block_id),
            "size": int(n),
            "internal_edges": int(e_rr),
            "external_edges": int(external_edges),
            "degree_in": int(k_in),
            "degree_out": int(k_out),
            "density": float(density),
            "q_r": float(q_r),
        }

    return block_stats


def attrib2graph(graph: gt.Graph, state: gt.NestedBlockState) -> gt.Graph:
    """
    Include the block attributes in the graph vertices as vertex properties.
    This function extracts block statistics from the SBM state and attaches them to the graph for further analysis or visualization.

    Args:
        graph (gt.Graph): The input graph-tool graph.
        state (gt.NestedBlockState): The nested block state from which to extract block statistics.
    Returns:
        gt.Graph: The graph with block attributes added as vertex properties.
    """
    state = state.copy()
    levels = state.get_levels()
    base_state = levels[0]  # Level 0 = finest partition
    blocks = base_state.get_blocks()

    block_array = blocks.a  # This is a numpy array

    # Creating the block property first
    unique_blocks = sorted(set(block_array))
    # B = len(unique_blocks) # Total number of blocks

    block_stats = get_block_attrib(
        graph=graph, unique_blocks=unique_blocks, block_array=block_array
    )

    # Define block properties
    block_id_prop = graph.new_vertex_property("int")
    block_size_prop = graph.new_vertex_property("int")
    block_internal_edges_prop = graph.new_vertex_property("int")
    block_external_edges_prop = graph.new_vertex_property("int")
    block_degree_in_prop = graph.new_vertex_property("int")
    block_degree_out_prop = graph.new_vertex_property("int")
    block_density_prop = graph.new_vertex_property("float")
    block_qr_prop = graph.new_vertex_property("float")

    # Include attributes in vertices
    for v in graph.vertices():
        block_id = block_array[int(v)]
        stats = block_stats[block_id]
        block_id_prop[v] = stats["block_id"]
        block_size_prop[v] = stats["size"]
        block_internal_edges_prop[v] = stats["internal_edges"]
        block_external_edges_prop[v] = stats["external_edges"]
        block_degree_in_prop[v] = stats["degree_in"]
        block_degree_out_prop[v] = stats["degree_out"]
        block_density_prop[v] = stats["density"]
        block_qr_prop[v] = stats["q_r"]

    # Assign properties to graph
    logger.info("📎 Attaching SBM metadata to the graph object…")

    graph.vertex_properties["SBM_l0_block_id"] = block_id_prop
    graph.vertex_properties["SBM_l0_block_size"] = block_size_prop
    graph.vertex_properties["SBM_l0_block_internal_edges"] = block_internal_edges_prop
    graph.vertex_properties["SBM_l0_block_external_edges"] = block_external_edges_prop
    graph.vertex_properties["SBM_l0_block_degree_in"] = block_degree_in_prop
    graph.vertex_properties["SBM_l0_block_degree_out"] = block_degree_out_prop
    graph.vertex_properties["SBM_l0_block_density"] = block_density_prop
    graph.vertex_properties["SBM_l0_block_qr"] = block_qr_prop

    logger.info("Block attributes added to graph vertices.")
    return graph


def export_graph(graph: gt.Graph, file_path: str):
    """
    Export the graph to a GraphML file.

    Args:
        graph (gt.Graph): The graph-tool graph to export.
        file_path (str): Path to the output GraphML file.
    """
    try:
        graph.save(file_path, fmt="graphml")
        logger.info(f"Graph exported successfully to {file_path}.")
    except Exception as e:
        logger.error(f"Failed to export graph to {file_path}: {e}")
        raise


def load_graph_tool_graph(edges_file: str) -> gt.Graph:
    # EDGES NETWORK FILE

    logger.info("Loading edges files ...")
    df = pd.read_csv(edges_file, sep="\t")
    logger.info("Edges loaded ...")

    g = gt.Graph(directed=False)

    # eprops: list of edge PropertyMaps to fill, matching column order after src/tgt
    e_pearson = g.new_edge_property("double")
    e_pval = g.new_edge_property("double")
    e_padj = g.new_edge_property("double")
    e_weight = g.new_edge_property("double")

    edges = df[["gene1", "gene2", "pearson", "pval", "padj", "weight"]].values

    logger.info("Adding edges into graph-tool element")
    vertex_names = g.add_edge_list(
        edges,
        hashed=True,  # gene1/gene2 are strings, not integers
        eprops=[e_pearson, e_pval, e_padj, e_weight],
    )

    g.vertex_properties["name"] = vertex_names
    g.edge_properties["pearson"] = e_pearson
    g.edge_properties["pval"] = e_pval
    g.edge_properties["padj"] = e_padj
    g.edge_properties["weight"] = e_weight
    
    logger.info("--- WEIGHTS ---")
    logger.info(f"PropertyMap type: {e_weight.value_type()}")
    
    # Acessa os dados brutos da propriedade como um array Numpy
    pesos_array = e_weight.a 
    
    logger.info(f"Lower: {np.min(pesos_array)}")
    logger.info(f"Upper: {np.max(pesos_array)}")
    logger.info(f"Are negatives? {np.any(pesos_array < 0)}")
    logger.info(f"Are zeros? {np.any(pesos_array == 0)}")
    logger.info(f"Are NaNs? {np.isnan(pesos_array).any()}")
    logger.info(f"Are Infs? {np.isinf(pesos_array).any()}")
    logger.info("-----------------------------")

    return g


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="""
SBM Inference CLI — Multi Layer Network
--------------------------------------

This tool performs Stochastic Block Model (SBM) inference on GraphML networks.
It supports this inference mode:

  • edge_analysis - Edge uncertainty analysis using weighted SBM

Typical workflow:
    1) Load a GraphML network using --input
    2) Tune annealing / MCMC parameters as needed
    3) Save output to a new GraphML using --output

Examples:
    # Basic default SBM inference
    python script.py -i edges.tsv -o sbm_out.graphml

    # Verbose MCMC debug output
    python script.py -i edges.tsv -o out.graphml -v

Notes:
    - For real data, recommended parameters are:
        --anneal-niter 100+
        --marginal-niter 100+
        --equilibrate-wait 100+
        --runs 20-100 (multi_anneal mode)
    - CPU-intensive steps benefit greatly from --cpus N
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=str,
        required=True,
        help="Path to input edges file.",
        metavar="FILE",
    )


    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Path to output GraphML file.",
        metavar="FILE",
    )

    parser.add_argument(
        "--runs",
        type=int,
        default=100,
        help="Number of independent annealing runs.",
    )

    parser.add_argument(
        "--anneal-niter",
        type=int,
        default=100,
        help="Number of annealing iterations per run (should be higher in production). 🔥",
    )

    parser.add_argument(
        "--equilibrate-wait",
        type=int,
        default=100,
        help="Steps to wait for MCMC equilibration. 🧘",
    )

    parser.add_argument(
        "--marginal-niter",
        type=int,
        default=100,
        help="Iterations for MCMC marginal sampling. 🎲",
    )

    parser.add_argument(
        "--mcmc-niter",
        type=int,
        default=10,
        help="Number of MCMC iterations for each step. ♻️",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=7,
        help="Random seed for reproducibility.",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        default=False,
        action="store_true",
        help="Verbose MCMC sampling output.",
    )

    return parser


def main():
    args = get_parser().parse_args()

    # Set seed for reproducibility
    np.random.seed(args.seed)
    try:
        gt.seed_rng(args.seed)  # for graph-tool randomization
    except Exception as e:
        logger.error(f"Error occurred while seeding random number generator: {e}")
        sys.exit(1)

    # Show all parameters in pretty table
    param_table = Table(show_header=True, header_style="bold magenta")
    param_table.add_column("Parameter", style="dim", width=30)
    param_table.add_column("Value", style="bold cyan")
    for arg in vars(args):
        param_table.add_row(arg, str(getattr(args, arg)))
    console.print(Panel(param_table, title="[bold magenta]SBM Inference Parameters[/]"))

    # Load graph
    graph = load_graph_tool_graph(args.input)

    state_weighted, _ = check_weighted_model(graph, graph.edge_properties["weight"])
    state = sbm_annealing(
        state_weighted,
        anneal_niter=args.anneal_niter,
        equilibrate_wait=args.equilibrate_wait,
        marginal_niter=args.marginal_niter,
        mcmc_niter=args.mcmc_niter,
        verbose=args.verbose,
    )

    # Edge uncertainty reconstruction
    marginal_graph = get_uncertainty_reconstruction(
        graph,
        initial_state=state,
        equilibrate_wait=args.equilibrate_wait,
        mcmc_niter=args.mcmc_niter,
        verbose=args.verbose,
    )

    # Attach marginal probabilities as edge property
    graph, df_results = integrate_reconstruction_results(
        graph,
        marginal_graph,
        block_state=state,
        prob_prop="eprob",
        prediction_threshold=0.001,
    )

    # Show block statistics
    # Add block attributes to graph
    graph = attrib2graph(graph, state)
    # Export graph with SBM attributes
    export_graph(graph, args.output)
    # Also export the DataFrame to CSV
    csv_output = args.output.replace(".graphml", "_edge_reconstruction_report.csv")
    df_results.to_csv(csv_output, index=False)
    logger.info(f"Edge reconstruction report exported to {csv_output}.")


if __name__ == "__main__":
    main()
