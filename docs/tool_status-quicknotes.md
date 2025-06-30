# Tool Status
_First test each module independently with example data from each tool's own repo_

- Lemur: working with example from repo
  - Need to include the optional parameters listed in `def parse_args` function [line 79](https://github.com/treangenlab/lemur/blob/main/lemur#L79)

- Magnet: Trying to make a conda env by packaging the named dependencies in readme (in the `dependencies.yml` file)
  - Micromamba makes an empty env for some reason ; troubleshoot this
- EMU: 
  - Copied example from EMU repo
  - Not working due to the `meta` input issue

Error:
```sh
Handling unexpected condition for
  task: name=EMU_ABUNDANCE; work-dir=null
  error [nextflow.exception.ProcessUnrecoverableException]: Not a valid path value type: groovyx.gpars.dataflow.DataflowBroadcast (DataflowBroadcast around DataflowStream[?])
```

# Process to make a nextflow module

1. Clone the tool's repo into a temp directory (outside this repo)
2. (optional) Test the tool with example data from the tool's own repo
3. Copy the module template from `modules/module_template.nf` to `modules/local/{tool_name}/main.nf`
4. Check for the tool's conda repo to call in the module-process's conda definition
5. For windsurf-AI's help in making the module, copy the tool's main script or readme of how to use it to the `modules/local/{tool_name}` directory as a placeholder file
  
