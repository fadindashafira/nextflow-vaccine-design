#!/bin/bash

# Run script for in silico vaccine design pipeline
# This script simplifies running the pipeline on a laptop

# Display usage information
show_usage() {
    echo "Usage: ./run_pipeline.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -a, --accession ACCESSION   NCBI protein accession number (default: BAL61230.1)"
    echo "  -o, --outdir DIRECTORY      Output directory (default: results)"
    echo "  -m, --run-md                Run molecular dynamics simulation (default: false)"
    echo "  -t, --test                  Run with test parameters"
    echo "  -g, --dag                   Generate DAG visualization"
    echo "  -r, --resume                Resume from the last successful run"
    echo "  -h, --help                  Show this help message"
    echo ""
    echo "Examples:"
    echo "  ./run_pipeline.sh"
    echo "  ./run_pipeline.sh -a 'QHD43416.1' -o 'covid_results'"
    echo "  ./run_pipeline.sh --accession 'BAL61230.1' --run-md --dag"
}

# Parse command line arguments
ACCESSION="BAL61230.1"
OUTDIR="results"
RUN_MD="false"
TEST_MODE="false"
GENERATE_DAG="true"
RESUME_MODE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -a|--accession)
            ACCESSION="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -m|--run-md)
            RUN_MD="true"
            shift
            ;;
        -t|--test)
            TEST_MODE="true"
            shift
            ;;
        -g|--dag)
            GENERATE_DAG="true"
            shift
            ;;
        -r|--resume)
            RESUME_MODE="true"
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            show_usage
            exit 1
            ;;
    esac
done

# Make the script executable
chmod +x "$0"

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "Error: Nextflow is not installed. Please install it first."
    echo "Run: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# Check if Graphviz is installed for DAG visualization
DAG_VISUALIZE="false"
if [ "$GENERATE_DAG" = "true" ] && command -v dag &> /dev/null; then
    DAG_VISUALIZE="true"
elif [ "$GENERATE_DAG" = "true" ]; then
    echo "Warning: Graphviz not installed. DAG visualization will be skipped."
fi

# Configure JVM memory for Nextflow
export NXF_OPTS="-Xms500M -Xmx2G"

# Start pipeline execution
echo "Starting in silico vaccine design pipeline..."
echo "Accession: $ACCESSION"
echo "Output directory: $OUTDIR"
echo "Run molecular dynamics: $RUN_MD"

# Prepare Nextflow command
if [ "$TEST_MODE" = "true" ]; then
    CMD="nextflow run main.nf -profile standard -with-report execution_report.html -with-timeline timeline.html"
elif [ "$RESUME_MODE" = "true" ]; then
    CMD="nextflow run main.nf -profile standard -resume"
else
    CMD="nextflow run main.nf -profile standard --accession '$ACCESSION' --outdir '$OUTDIR' --run_md $RUN_MD"
fi

# Add DAG generation to command if requested
if [ "$GENERATE_DAG" = "true" ]; then
    CMD+=" -with-dag $OUTDIR/pipeline_dag.dag"
fi

# Execute the pipeline
echo "Executing: $CMD"
eval $CMD

# Store the execution status
EXEC_STATUS=$?

# Generate DAG visualizations if requested and successful
if [ "$EXEC_STATUS" -eq 0 ] && [ "$DAG_VISUALIZE" = "true" ]; then
    echo "Generating DAG visualizations..."
    
    # Ensure output directory exists
    mkdir -p "$OUTDIR"
    
    # Convert to PDF
    dot -Tpdf "$OUTDIR/pipeline_dag.dag" -o "$OUTDIR/pipeline_dag.pdf"
    
    # Convert to PNG
    dot -Tpng "$OUTDIR/pipeline_dag.dag" -o "$OUTDIR/pipeline_dag.png"
    
    echo "DAG visualizations created:"
    echo "- $OUTDIR/pipeline_dag.pdf"
    echo "- $OUTDIR/pipeline_dag.png"
fi

# Check execution status
if [ $EXEC_STATUS -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
    echo "Results are available in: $OUTDIR"
    echo ""
    echo "Key output files:"
    echo "- Vaccine construct: $OUTDIR/vaccine_construct.fasta"
    echo "- Evaluation report: $OUTDIR/evaluation/vaccine_evaluation.txt"
    
    if [ "$RUN_MD" = "true" ]; then
        echo "- MD simulation report: $OUTDIR/molecular_dynamics/md_report.txt"
    fi
else
    echo ""
    echo "Pipeline execution failed. Check the logs for details."
    echo "To resume from the last successful step, run: ./run_pipeline.sh -r"
fi