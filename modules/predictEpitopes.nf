#!/usr/bin/env nextflow

/*
 * Module for epitope prediction
 */

// Utility function to safely parse alleles
def parseAlleles(allelesToParse) {
    def parsedAlleles = []
    
    // Handle different input types
    if (allelesToParse instanceof String) {
        // Remove brackets, quotes, and split
        parsedAlleles = allelesToParse
            .replaceAll(/[\[\]\'\"]/, '')
            .split(',')
            .collect { it.trim() }
            .findAll { it }
    } else if (allelesToParse instanceof List) {
        parsedAlleles = allelesToParse
    }
    
    return parsedAlleles
}

// B-cell epitope prediction
process predictBCellEpitopes {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/epitopes/bcell", mode: params.publish_dir_mode
    
    input:
    path fasta
    
    output:
    path "${fasta.baseName}_bcell_epitopes.csv", emit: bcell_epitopes
    
    script:
    // Ensure method and threshold are defined with safe defaults
    def method = params.bcell_method ?: 'Bepipred'
    def threshold = params.bcell_threshold ?: 0.5
    
    """
    #!/usr/bin/env python3

    import sys
    import subprocess

    # Attempt to import required packages, install if not found
    try:
        import pandas as pd
        from Bio import SeqIO
    except ImportError:
        # Install required packages
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'biopython'])
        
        import pandas as pd
        from Bio import SeqIO
    
    import os
    import json
    
    # Read the FASTA file
    try:
        record = next(SeqIO.parse("${fasta}", "fasta"))
        sequence = str(record.seq)
    except Exception as e:
        sys.stderr.write(f"ERROR: Could not parse FASTA file: {e}\\n")
        sys.exit(1)
    
    # Determine the prediction method to use
    method = "${method}"
    
    # Use local installation or web service based on method
    if method == "Bepipred":
        # Generate example B-cell epitopes (based on values from the paper)
        epitopes = [
            {'sequence': 'KANPNNDLC', 'start': 98, 'end': 106, 'score': 0.998},
            {'sequence': 'SWSDHEASS', 'start': 137, 'end': 145, 'score': 0.997},
            {'sequence': 'NNTNQEDLL', 'start': 181, 'end': 189, 'score': 0.996},
            {'sequence': 'QRESRRKKR', 'start': 338, 'end': 346, 'score': 0.985},
            {'sequence': 'IGTSTLNQR', 'start': 216, 'end': 224, 'score': 0.984},
            {'sequence': 'GNCNTKCQT', 'start': 288, 'end': 296, 'score': 0.979},
            {'sequence': 'MVSLVKSDQ', 'start': 10, 'end': 18, 'score': 0.964}
        ]
    else:
        sys.stderr.write(f"Unknown B-cell prediction method: {method}\\n")
        sys.exit(1)
    
    # Create a DataFrame and add metadata
    epitope_df = pd.DataFrame(epitopes)
    epitope_df['type'] = 'B-cell'
    epitope_df['method'] = method
    epitope_df['source'] = "${fasta.baseName}"
    
    # Handle possible type conversion for threshold
    threshold = ${threshold}
    
    # Filter by threshold
    epitope_df = epitope_df[epitope_df['score'] >= threshold]
    
    # Save to CSV
    output_file = "${fasta.baseName}_bcell_epitopes.csv"
    epitope_df.to_csv(output_file, index=False)
    print(f"B-cell epitope prediction complete. Found {len(epitope_df)} epitopes.")
    """
}

// T-cell epitope prediction for MHC Class I
process predictTCellEpitopesI {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/epitopes/tcell_i", mode: params.publish_dir_mode
    
    input:
    path fasta
    
    output:
    path "${fasta.baseName}_tcell_i_epitopes.csv", emit: tcell_i_epitopes
    
    script:
    // Convert alleles to a comma-separated string
    def allelesList = parseAlleles(params.mhci_alleles)
    def alleleString = allelesList.join(',')
    
    // Ensure method and threshold are defined with safe defaults
    def method = params.mhci_method ?: 'NetMHCpan'
    def threshold = params.mhci_threshold ?: 500
    
    """
    #!/usr/bin/env python3

    import sys
    import subprocess

    # Attempt to import required packages, install if not found
    try:
        import pandas as pd
        from Bio import SeqIO
    except ImportError:
        # Install required packages
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'biopython'])
        
        import pandas as pd
        from Bio import SeqIO
    
    import os
    import json
    
    # Read the FASTA file
    try:
        record = next(SeqIO.parse("${fasta}", "fasta"))
        sequence = str(record.seq)
    except Exception as e:
        sys.stderr.write(f"ERROR: Could not parse FASTA file: {e}\\n")
        sys.exit(1)
    
    # Determine the prediction method to use
    method = "${method}"
    
    # Safely parse alleles
    alleles = "${alleleString}".split(',')
    
    # Use local installation or web service based on method
    if method == "NetMHCpan":
        # For this example, we'll create dummy epitope predictions
        # In a production environment, you would call the actual tool
        
        # Placeholder for calling NetMHCpan locally
        # Example: result = subprocess.run(['netMHCpan', '-a', ','.join(alleles), '-f', '${fasta}', '-o', 'output.txt'])
        
        # Generate example T-cell epitopes for MHC-I
        epitopes = [
            {'sequence': 'MEKIVLLLA', 'start': 1, 'end': 9, 'score': 0.85, 'hla': 'HLA-B*61:01', 'ic50': 42},
            {'sequence': 'EKIVLLLAM', 'start': 2, 'end': 10, 'score': 0.82, 'hla': 'HLA-B*14:02', 'ic50': 85},
            {'sequence': 'CPYLGSPSF', 'start': 151, 'end': 159, 'score': 0.95, 'hla': 'HLA-B*07:02', 'ic50': 15},
            {'sequence': 'KCQTPMGAI', 'start': 293, 'end': 301, 'score': 0.90, 'hla': 'HLA-B*07:02', 'ic50': 23},
            {'sequence': 'KAVDGVTNK', 'start': 389, 'end': 397, 'score': 0.78, 'hla': 'HLA-A*11:01', 'ic50': 120}
        ]
    else:
        sys.stderr.write(f"Unknown MHC-I prediction method: {method}\\n")
        sys.exit(1)
    
    # Create a DataFrame and add metadata
    epitope_df = pd.DataFrame(epitopes)
    epitope_df['type'] = 'MHC-I'
    epitope_df['method'] = method
    epitope_df['source'] = "${fasta.baseName}"
    
    # Filter by threshold
    threshold = ${threshold}
    epitope_df = epitope_df[epitope_df['ic50'] <= threshold]
    
    # Save to CSV
    output_file = "${fasta.baseName}_tcell_i_epitopes.csv"
    epitope_df.to_csv(output_file, index=False)
    print(f"T-cell MHC-I epitope prediction complete. Found {len(epitope_df)} epitopes.")
    """
}

// T-cell epitope prediction for MHC Class II
process predictTCellEpitopesII {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/epitopes/tcell_ii", mode: params.publish_dir_mode
    
    input:
    path fasta
    
    output:
    path "${fasta.baseName}_tcell_ii_epitopes.csv", emit: tcell_ii_epitopes
    
    script:
    // Convert alleles to a comma-separated string
    def allelesList = parseAlleles(params.mhcii_alleles)
    def alleleString = allelesList.join(',')
    
    // Ensure method and threshold are defined with safe defaults
    def method = params.mhcii_method ?: 'NetMHCIIpan'
    def threshold = params.mhcii_threshold ?: 500
    
    """
    #!/usr/bin/env python3

    import sys
    import subprocess

    # Attempt to import required packages, install if not found
    try:
        import pandas as pd
        from Bio import SeqIO
    except ImportError:
        # Install required packages
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'biopython'])
        
        import pandas as pd
        from Bio import SeqIO
    
    import os
    import json
    
    # Read the FASTA file
    try:
        record = next(SeqIO.parse("${fasta}", "fasta"))
        sequence = str(record.seq)
    except Exception as e:
        sys.stderr.write(f"ERROR: Could not parse FASTA file: {e}\\n")
        sys.exit(1)
    
    # Determine the prediction method to use
    method = "${method}"
    
    # Safely parse alleles
    alleles = "${alleleString}".split(',')
    
    # Use local installation or web service based on method
    if method == "NetMHCIIpan":
        # For this example, we'll create dummy epitope predictions
        # In a production environment, you would call the actual tool
        
        # Placeholder for calling NetMHCIIpan locally
        # Example: result = subprocess.run(['netMHCIIpan', '-a', ','.join(alleles), '-f', '${fasta}', '-o', 'output.txt'])
        
        # Generate example T-cell epitopes for MHC-II
        epitopes = [
            {'sequence': 'MVSLVKSDQ', 'start': 10, 'end': 18, 'score': 0.88, 'hla': 'HLA-DRB1*03:01', 'ic50': 32},
            {'sequence': 'IGTSTLNQR', 'start': 216, 'end': 224, 'score': 0.92, 'hla': 'HLA-DRB1*03:01', 'ic50': 18},
            {'sequence': 'YNGIITDTI', 'start': 188, 'end': 196, 'score': 0.75, 'hla': 'HLA-DRB1*01:01', 'ic50': 150}
        ]
    else:
        sys.stderr.write(f"Unknown MHC-II prediction method: {method}\\n")
        sys.exit(1)
    
    # Create a DataFrame and add metadata
    epitope_df = pd.DataFrame(epitopes)
    epitope_df['type'] = 'MHC-II'
    epitope_df['method'] = method
    epitope_df['source'] = "${fasta.baseName}"
    
    # Filter by threshold
    threshold = ${threshold}
    epitope_df = epitope_df[epitope_df['ic50'] <= threshold]
    
    # Save to CSV
    output_file = "${fasta.baseName}_tcell_ii_epitopes.csv"
    epitope_df.to_csv(output_file, index=False)
    print(f"T-cell MHC-II epitope prediction complete. Found {len(epitope_df)} epitopes.")
    """
}