# Interactive HTML Reporting Framework for Somatem

## Overview
This document outlines the framework for implementing interactive HTML reporting in the Somatem pipeline, leveraging Nextflow's built-in reporting capabilities and MultiQC integration for a comprehensive, interactive pipeline report. The implementation ensures both live monitoring during pipeline execution and full functionality when downloaded for offline/archival use.

## Report Persistence
The reporting system generates self-contained HTML reports that work both during pipeline execution and when downloaded for offline viewing:

### Live Execution Features
- Real-time process status updates via weblog
- Live DAG updates showing current execution state
- Progressive MultiQC report generation

### Offline/Archived Features
- Complete DAG visualization with all process states
- Full MultiQC report with all tool outputs
- Interactive data exploration capabilities
- All JavaScript and styling bundled in HTML
- No external dependencies required for viewing

## Components

### 1. Core Technologies
- **Nextflow Built-in Reports** (`-with-dag`, `-with-report`, `-with-timeline`)
- **MultiQC** - For aggregating and visualizing tool outputs
- **JavaScript/jQuery** - For enhanced interactivity
- **Mermaid.js** - For additional workflow visualization options
- **Server-Sent Events (SSE)** - For live updates via `-with-weblog`

### 2. Implementation Structure

```plaintext
assets/
└── report_templates/
    ├── custom_multiqc_config.yaml  # MultiQC configuration
    ├── dag_enhancer.js            # DAG interactivity scripts
    ├── custom_report.html         # Custom report template
    └── styles/
        └── custom_report.css      # Custom styling
```

### 3. Key Features

#### 3.1 Enhanced DAG Visualization
- Utilize Nextflow's built-in DAG visualization (`-with-dag`)
- Enhance DAG with custom JavaScript for node clicking
- Link DAG nodes to corresponding MultiQC sections
- Color-coded process states using Nextflow's execution data

```javascript
// Example flowchart structure in flowchart.js
{
  "nodes": [
    {"id": "input", "label": "Input FastQ"},
    {"id": "qc", "label": "NanoPlot QC"},
    {"id": "hostile", "label": "Host Removal"},
    // ... other processes
  ],
  "links": [
    {"source": "input", "target": "qc"},
    {"source": "qc", "target": "hostile"},
    // ... process connections
  ]
}
```

#### 3.2 MultiQC Integration
- Comprehensive tool output aggregation
- Built-in support for many bioinformatics tools
- Customizable report sections
- Interactive data visualization
- Support for custom content and statistics

#### 3.3 Live Updates Implementation

```javascript
// Using Server-Sent Events with Nextflow's weblog
const evtSource = new EventSource('http://localhost:PORT/events');

evtSource.onmessage = (event) => {
    const data = JSON.parse(event.data);
    switch(data.event) {
        case 'process_submitted':
            updateDagNode(data.process, 'running');
            break;
        case 'process_completed':
            updateDagNode(data.process, 'complete');
            updateMultiQCSection(data.process);
            break;
        case 'error':
            updateDagNode(data.process, 'failed');
            showError(data.error);
            break;
    }
};

// Example of DAG node enhancement
function enhanceDagNode(nodeId, outputs) {
    const node = document.querySelector(`#dag-${nodeId}`);
    node.addEventListener('click', () => {
        // Scroll to corresponding MultiQC section
        document.querySelector(`#multiqc-${nodeId}`).scrollIntoView();
        // Show process-specific outputs
        showOutputs(outputs);
    });
}
```

### 4. Integration with Nextflow

#### 4.1 Process Output Configuration
```nextflow
process EXAMPLE_PROCESS {
    publishDir "${params.outdir}/process_name", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.html')) "qc_reports/$filename"
            else if (filename.endsWith('.txt')) "stats/$filename"
            else if (filename.endsWith('.log')) "logs/$filename"
            else filename
        }
    
    script:
    """
    # Process commands here
    # Ensure output formats are MultiQC-compatible
    """
}
```

#### 4.2 MultiQC Process
```nextflow
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'  // Collect all QC reports
    path config from ch_multiqc_config

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . \
        --config $config \
        --interactive \
        -o . \
        --title "${params.title}" \
        --comment "Somatem pipeline report"
    """
}
```

### 5. Implementation Steps

1. **Configure Nextflow Reporting**
   ```bash
   nextflow run main.nf \
       -with-dag flowchart.html \
       -with-report execution_report.html \
       -with-timeline timeline.html \
       -with-weblog 'http://localhost:8000/events'
   ```

2. **Setup MultiQC Configuration**
   - Create custom MultiQC config
   - Define module order
   - Configure custom content
   - Set up custom visualization options

3. **Enhance DAG Visualization**
   - Add custom JavaScript for node interactions
   - Link nodes to MultiQC sections
   - Implement live updates using weblog events

4. **Process Output Integration**
   - Ensure tool outputs are MultiQC-compatible
   - Structure output directories for proper collection
   - Add custom MultiQC modules if needed

5. **Testing**
   - Validate DAG interactivity
   - Check MultiQC report generation
   - Test live update functionality
   - Verify mobile compatibility

### 6. Configuration

Add to `nextflow.config`:
```nextflow
params {
    // MultiQC settings
    multiqc_config = "$projectDir/assets/report_templates/custom_multiqc_config.yaml"
    max_multiqc_email_size = '25.MB'
    
    // Report customization
    title = 'Somatem Pipeline Report'
    custom_logo = "$projectDir/assets/logo.png"
    
    // Web report settings
    weblog_port = 8000
    enable_dag_visualization = true
    
    // Offline report settings
    bundle_resources = true     // Bundle all resources into self-contained HTML
    save_intermediate = true    // Save intermediate files for offline browsing
}

// Configure reporting
dag {
    enabled = true
    file = "${params.outdir}/pipeline_dag.html"
    overwrite = true
}

report {
    enabled = true
    file = "${params.outdir}/execution_report.html"
}

timeline {
    enabled = true
    file = "${params.outdir}/execution_timeline.html"
}

weblog {
    enabled = true
    url = "http://localhost:${params.weblog_port}/events"
}
```

### 7. Usage

Run the pipeline with reporting enabled:

```bash
nextflow run main.nf \
    -param-file assets/custom_metadata.yaml \
    -with-dag \
    -with-report \
    -with-timeline \
    -with-weblog
```

Reports will be available at:
- Pipeline DAG: `${params.outdir}/pipeline_dag.html`
- Execution Report: `${params.outdir}/execution_report.html`
- Timeline: `${params.outdir}/execution_timeline.html`
- MultiQC Report: `${params.outdir}/multiqc/multiqc_report.html`

## Security Considerations

- Use local-only access for weblog server
- Implement proper file permissions for output directories
- Consider access controls for sensitive data
- Use secure protocols if exposing reports externally

## Dependencies

Add to `nf_base_env.yml`:
```yaml
dependencies:
  - python>=3.8
  - multiqc>=1.14
  - jinja2>=3.0
  - markdown>=3.3
  - pymdown-extensions>=9.4
```