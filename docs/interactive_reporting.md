# Interactive HTML Reporting Framework for Somatem

## Overview
This document outlines the framework for implementing interactive HTML reporting in the Somatem pipeline, leveraging Nextflow's built-in reporting capabilities and MultiQC integration for a comprehensive, interactive pipeline report.

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
    ├── index.html          # Main dashboard template
    ├── flowchart.js        # D3.js pipeline visualization
    ├── live_updates.js     # WebSocket handling
    └── styles/
        └── dashboard.css   # Styling for the dashboard
```

### 3. Key Features

#### 3.1 Live Pipeline Flowchart
- Interactive D3.js visualization
- Color-coded process states (running, completed, failed)
- Clickable nodes linking to process details
- Auto-updates via WebSocket connection

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

#### 3.2 Process Output Tabs
- Dynamic tab generation for each process
- Automatic content type detection
- Support for:
  - Tables (CSV, TSV)
  - Plots (PNG, PDF, SVG)
  - Interactive visualizations (HTML)
  - Log files
  - MultiQC reports

#### 3.3 Live Updates Implementation

```javascript
// WebSocket connection for live updates
const ws = new WebSocket('ws://localhost:PORT');

ws.onmessage = (event) => {
  const update = JSON.parse(event.data);
  switch(update.type) {
    case 'process_start':
      updateFlowchart(update);
      break;
    case 'process_complete':
      updateProcessTab(update);
      break;
    case 'new_output':
      addOutputToTab(update);
      break;
  }
};
```

### 4. Integration with Nextflow

#### 4.1 Process Output Publishing
```nextflow
process EXAMPLE_PROCESS {
    publishDir "${params.outdir}/process_name", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.html')) "reports/$filename"
            else if (filename.endsWith('.csv')) "tables/$filename"
            else if (filename.endsWith('.png')) "plots/$filename"
            else filename
        }
    
    script:
    """
    # Process commands here
    """
}
```

#### 4.2 Report Generation Process
```nextflow
process GENERATE_REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    // Collect all process outputs

    output:
    path "interactive_report.html"

    script:
    """
    # Generate final HTML report
    # Combine all outputs
    # Setup WebSocket server
    """
}
```

### 5. Implementation Steps

1. **Setup Report Server**
   - Create a lightweight server (Node.js/Python) to handle WebSocket connections
   - Implement file watching for output directories
   - Setup route handling for different output types

2. **Create Base Templates**
   - Design main dashboard layout
   - Implement tab navigation system
   - Create flowchart visualization

3. **Process Integration**
   - Modify process definitions to publish outputs in structured format
   - Implement output processors for different file types
   - Add metadata generation for each process

4. **Testing**
   - Verify real-time updates
   - Test different output types
   - Validate mobile responsiveness

### 6. Configuration

Add to `nextflow.config`:
```nextflow
params {
    // Interactive report settings
    interactive_report = true
    report_port = 8080
    report_update_interval = 10 // seconds
    report_dir = "${params.outdir}/interactive_report"
}
```

### 7. Usage

To enable interactive reporting:

```bash
nextflow run main.nf -param-file assets/custom_metadata.yaml --interactive_report true
```

The interactive report will be available at:
```
http://localhost:8080
```

## Security Considerations

- Implement authentication for sensitive data
- Use HTTPS for secure connections
- Sanitize all user inputs
- Implement rate limiting for WebSocket connections
- Configure proper CORS policies

## Dependencies

Add to `nf_base_env.yml`:
```yaml
dependencies:
  - nodejs>=14.0.0
  - python>=3.8
  - websockets
  - aiohttp
  - plotly
```