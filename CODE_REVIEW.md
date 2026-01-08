# Somatem Pipeline - Code Review Report

**Date:** 2026-01-08  
**Reviewer:** Seqera AI  
**Pipeline Version:** 1.0.0dev  
**Branch:** global_publish_i41

---

## Executive Summary

This code review evaluates the Somatem long-read metagenomics pipeline for code quality, Nextflow best practices, and opportunities for improvement. The pipeline is well-structured and demonstrates good modular design. Key recommendations focus on:

1. **Migrating to global publishing syntax** (v25.10+ output blocks)
2. **Standardizing module patterns** 
3. **Improving error handling and validation**
4. **Enhancing documentation**
5. **Optimizing resource allocation**

**Overall Assessment:** ⭐⭐⭐⭐ (4/5) - Strong foundation with room for modernization

---

## 1. Publishing Strategy - PRIORITY REFACTOR

### Current State
- **Mixed approach:** Combination of process-level `publishDir` and workflow-level `output` blocks
- **8 modules** have hardcoded `publishDir` directives
- **Basic global outputs** in main.nf (mapping, binning_tables, binning_fasta)

### Issues
❌ **Inconsistent publishing patterns** across modules  
❌ **Hardcoded paths** make customization difficult  
❌ **Limited use of v25.10+ output blocks**  
❌ **Difficult to track what gets published** where  

### Recommendations ✅

**HIGH PRIORITY:** Migrate to global publishing using Nextflow v25.10+ `output` blocks

#### Benefits:
- Centralized control of all pipeline outputs
- Clear visibility of what gets published
- Dynamic path configuration using closures
- Easier to maintain and modify
- Better separation of concerns (modules handle processing, main workflow handles publishing)

#### Implementation Strategy:

1. **Remove all `publishDir` directives from modules:**
   - `modules/nf-core/flye/main.nf`
   - `modules/local/bakta/bakta/main.nf`
   - `modules/local/singlem/pipe/main.nf`
   - `modules/local/singlem/appraise/main.nf`
   - `modules/local/checkm2/parse/main.nf`
   - `modules/local/taxburst/main.nf`
   - `modules/nf-core/checkm2/predict/main.nf`
   - `modules/nf-core/samtools/coverage/main.nf`

2. **Expand output block in main.nf** to include all key outputs

3. **Use descriptive path configurations** with sample-aware organization

---

## 2. Code Structure & Organization

### Strengths ✅
- **Excellent modular design** with clear separation of subworkflows
- **Logical workflow organization:** preprocessing → profiling/assembly → analysis
- **Good use of channels** and data flow patterns
- **Comprehensive metadata handling** with meta maps

### Areas for Improvement

#### 2.1 Channel Operations

**Issue:** Complex channel transformations in `assembly_mags.nf`
```groovy
// Lines 207-225: Complex metadata manipulation
ch_bins_for_annotation = SEMIBIN_SINGLEEASYBIN.out.output_fasta
    .transpose()
    .map { meta, bin ->
        def new_meta = meta.clone()
        def bin_filename = bin.name
        def clean_bin_name = bin_filename.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
        new_meta.id = clean_bin_name
        new_meta.sample_id = meta.id
        // ... more transformations
    }
```

**Recommendation:**
- Extract channel transformation logic into **helper functions**
- Add **inline documentation** explaining transformations
- Consider creating a **utils module** for common operations

#### 2.2 Conditional Logic

**Issue:** Analysis-type conditionals scattered across workflows
```groovy
// In somatem.nf
if (params.analysis_type == "taxonomic-profiling") { ... }
if (params.analysis_type == "assembly") { ... }
if (params.analysis_type == "genome-dynamics") { ... }
```

**Recommendation:**
- Consider using **sub-workflows as optional includes**
- Use **workflow handlers** for cleaner conditional execution
- Validate `analysis_type` parameter early in pipeline initialization

---

## 3. Error Handling & Validation

### Current State
✅ Good validation in BAKTA_BAKTA module (FASTA checks)  
✅ Flye mode validation  
⚠️ Limited parameter validation elsewhere  

### Recommendations

#### 3.1 Input Validation
```groovy
// Add to main workflow or PIPELINE_INITIALISATION
workflow {
    // Validate analysis_type parameter
    def valid_analysis_types = ['taxonomic-profiling', 'assembly', 'genome-dynamics']
    if (!valid_analysis_types.contains(params.analysis_type)) {
        error "Invalid analysis_type: ${params.analysis_type}. Must be one of: ${valid_analysis_types.join(', ')}"
    }
    
    // Validate data_type
    def valid_data_types = ['metagenomics', '16S']
    if (!valid_data_types.contains(params.data_type)) {
        error "Invalid data_type: ${params.data_type}. Must be one of: ${valid_data_types.join(', ')}"
    }
}
```

#### 3.2 Database Validation
Add checks for required databases before processing:
```groovy
// In DOWNLOAD_DBS or at workflow start
if (params.analysis_type == "assembly" && !file(params.checkm2_db).exists()) {
    log.warn "CheckM2 database not found. Will download..."
}
```

#### 3.3 Empty Channel Handling
Add protective checks for empty channels:
```groovy
// Example pattern
PREPROCESSING.out.clean_reads
    .ifEmpty { error "No reads passed quality filtering" }
    | NEXT_PROCESS
```

---

## 4. Nextflow DSL2 Best Practices

### 4.1 Modern Syntax (v25.10+)

**Current Issues:**
❌ Comments suggest using `Channel.from()` which is deprecated  
❌ Not using `channel.of()` namespace consistently  

**Recommendations:**
```groovy
// OLD (deprecated)
Channel.from('Hello', 'World')

// NEW (v25.10+)
channel.of('Hello', 'World')
```

### 4.2 Workflow Output Handlers

**Recommendation:** Use modern workflow handlers in entry workflow
```groovy
workflow {
    main:
    // workflow logic
    
    publish:
    // outputs
    
    onComplete:
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'FAILED'}"
    log.info "Results published to: ${params.outdir}"
    
    onError:
    log.error "Pipeline execution failed"
    log.error "Error message: ${workflow.errorMessage}"
}
```

### 4.3 Process Directives

**Issue:** Mixed use of labels and explicit resource allocation
```groovy
// In modules.config
withName: 'HOSTILE_CLEAN' {
    cpus = 3
    memory = '18.GB' 
}
```

**Recommendation:** Prefer labels and only override when necessary
```groovy
withName: 'HOSTILE_CLEAN' {
    label = 'process_medium_memory'  // Use existing or create new label
}
```

---

## 5. Configuration Management

### Strengths ✅
- Comprehensive parameter definitions
- Good use of profiles
- Clear separation of module-specific configs

### Issues & Recommendations

#### 5.1 Commented Code
**Issue:** Large blocks of commented resource allocations (lines 176-245 in modules.config)

**Recommendation:**
- Remove obsolete commented code or move to separate documentation
- If keeping for reference, add clear explanation of why it's commented

#### 5.2 Database Paths
**Issue:** Hardcoded database base directory
```groovy
params.db_base_dir = "/home/dbs"
```

**Recommendation:**
- Use environment variables or automatic detection
- Provide clear documentation for setting up databases
```groovy
params.db_base_dir = System.getenv('SOMATEM_DB_DIR') ?: "${projectDir}/databases"
```

#### 5.3 Parameter Organization
**Recommendation:** Group related parameters more clearly
```groovy
// Preprocessing parameters
params {
    preprocessing {
        maxlength = 30000
        minlen = 250
        minq = 10
    }
}
```

---

## 6. Module-Specific Issues

### 6.1 FLYE Module
**File:** `modules/nf-core/flye/main.nf`

**Issues:**
- ❌ Hardcoded `publishDir` directives (lines 6-8)
- ⚠️ Inconsistent compression (some outputs gzipped, log not)

**Recommendations:**
- Remove `publishDir` and use global output block
- Standardize compression approach
- Consider making log compression optional for easier debugging

### 6.2 BAKTA_BAKTA Module
**File:** `modules/local/bakta/bakta/main.nf`

**Strengths:** ✅ Excellent error handling and validation

**Issues:**
- ❌ Hardcoded `publishDir` with nested sample structure
- ⚠️ Creates empty files on error (may hide failures)

**Recommendations:**
```groovy
// Instead of creating empty files, consider:
if [ ! -f "\$file" ]; then
    echo "ERROR: Required output \$file not generated"
    exit 1  // Fail fast rather than creating empty files
fi
```

### 6.3 SINGLEM_PIPE Module
**File:** `modules/local/singlem/pipe/main.nf`

**Issues:**
- ❌ Hardcoded `publishDir`
- ⚠️ Complex conditional output handling
- ⚠️ Potential for inconsistent outputs based on conditions

**Recommendations:**
- Simplify output handling logic
- Use consistent output naming patterns
- Document expected outputs for each sample type

### 6.4 Channel Joining Complexity
**File:** `subworkflows/local/assembly_mags.nf` (lines 238-252)

**Issue:** Complex completeness map joining logic
```groovy
ch_bins_with_completeness = ch_bins_for_annotation
    .map { meta, bin -> 
        def join_key = meta.sample_id
        [join_key, meta, bin]
    }
    .combine(
        ch_completeness_with_meta.map { meta, completeness_map -> 
            [meta.id, completeness_map] 
        }, 
        by: 0
    )
```

**Recommendation:**
- Extract to a dedicated function
- Add comprehensive comments
- Consider using `.join()` operator instead of `.combine()` where appropriate

---

## 7. Documentation

### Current State
✅ Good inline comments in complex sections  
✅ Process tags for tracking  
✅ Log messages for key steps  
⚠️ Limited high-level documentation  

### Recommendations

#### 7.1 Module Documentation
Add standard documentation blocks to all modules:
```groovy
/*
 * Module: PROCESS_NAME
 * Description: Brief description of what this process does
 * Input:
 *   - meta: Map containing sample metadata
 *   - reads: Path to input files
 * Output:
 *   - results: Description of output
 * Notes:
 *   - Any special considerations
 */
process PROCESS_NAME {
    ...
}
```

#### 7.2 Workflow Documentation
Create or enhance:
- **Usage documentation** for each analysis type
- **Parameter descriptions** with valid values
- **Output structure** documentation
- **Troubleshooting guide**

#### 7.3 Code Comments
Add explanatory comments for:
- Complex channel transformations
- Non-obvious logic decisions
- Workarounds for known issues

---

## 8. Testing & Quality Assurance

### Current State
✅ Test configurations (test.config, test_full.config)  
✅ Stub sections in modules  
⚠️ Limited visibility into test coverage  

### Recommendations

#### 8.1 Add nf-test
```groovy
// Example test structure
nextflow_process {
    name "Test FLYE assembly"
    script "modules/nf-core/flye/main.nf"
    process "FLYE"
    
    test("Should assemble reads") {
        when {
            params {
                flye_mode = "nano-hq"
            }
            process {
                """
                input[0] = [ [id: 'test'], file('test.fastq.gz') ]
                input[1] = 'nano-hq'
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.fasta).exists()
        }
    }
}
```

#### 8.2 Continuous Integration
- Add GitHub Actions workflow for automated testing
- Run linting on pull requests
- Test multiple profiles (docker, singularity, conda)

#### 8.3 Validation Checks
Add output validation:
- Check file formats (FASTA, BAM, etc.)
- Validate completeness thresholds are met
- Verify required outputs exist

---

## 9. Performance Optimization

### 9.1 Resource Allocation

**Current Issues:**
- Some processes override label resources (may be intentional)
- No dynamic resource allocation based on input size

**Recommendations:**
```groovy
// Dynamic memory allocation based on input
withName: 'FLYE' {
    memory = { 36.GB * task.attempt * (meta.read_count ? Math.ceil(meta.read_count / 1000000) : 1) }
}
```

### 9.2 Parallelization

**Observations:**
✅ Good use of `maxForks` for problematic processes (FinalNanoPlot)  
⚠️ Some processes might benefit from increased parallelization  

**Recommendations:**
- Profile pipeline execution to identify bottlenecks
- Consider splitting large processes into parallel chunks where possible
- Use `executor.queueSize` to optimize job submission

### 9.3 Caching

**Good practices observed:**
✅ Using `storeDir` for database downloads  
✅ Proper work directory management  

**Recommendations:**
- Document which processes are cached
- Add cache strategy directives where appropriate:
```groovy
cache = 'lenient'  // For processes that may have filesystem timestamp issues
```

---

## 10. Security & Reproducibility

### 10.1 Container Usage
**Current:** Multiple container systems supported (Docker, Singularity, Apptainer)

**Recommendations:**
✅ Good support for multiple systems  
➕ Consider Wave for automatic container building  
➕ Use container digests for exact reproducibility:
```groovy
container = 'quay.io/biocontainers/flye:2.9.6--py311h6a68c12_0@sha256:abc123...'
```

### 10.2 Version Tracking
✅ Good version emission from processes  
➕ Consider adding git commit hash to output metadata  
➕ Generate comprehensive methods description  

---

## 11. Priority Action Items

### Immediate (High Priority)
1. ✅ **Implement global publishing syntax** (addresses main request)
2. ✅ **Remove hardcoded publishDir from modules**
3. ✅ **Add input parameter validation**
4. ✅ **Standardize output compression**

### Short Term (Medium Priority)
5. Clean up commented code in modules.config
6. Extract complex channel operations to helper functions
7. Add workflow completion handlers
8. Improve error messages for common failures

### Long Term (Lower Priority)
9. Implement comprehensive nf-test suite
10. Add performance profiling and optimization
11. Create detailed user documentation
12. Set up CI/CD pipeline

---

## 12. Positive Highlights ⭐

The following aspects of the pipeline demonstrate excellent practice:

1. **Modular Design:** Clean separation of concerns across subworkflows
2. **Metadata Handling:** Comprehensive use of meta maps for sample tracking
3. **Conditional Processing:** Smart handling of different analysis types
4. **Completeness Filtering:** Intelligent filtering of bins based on quality metrics
5. **Database Management:** Automated database download and caching
6. **Error Handling in BAKTA:** Exemplary validation and error checking
7. **Resource Management:** Thoughtful allocation with retry multipliers
8. **Logging:** Helpful progress messages throughout pipeline
9. **Profile Support:** Comprehensive support for different execution environments
10. **Scientific Rigor:** Appropriate tool selection and parameterization

---

## 13. Example Refactored Code

### Before (Current):
```groovy
// In FLYE module
publishDir "${params.output_dir}/assembly/${meta.id}", mode: 'copy', pattern: "*.fasta"
publishDir "${params.output_dir}/assembly/${meta.id}", mode: 'copy', pattern: "*.gfa"
publishDir "${params.output_dir}/assembly/${meta.id}", mode: 'copy', pattern: "*.log"
```

### After (Recommended):
```groovy
// In main.nf
output {
    assembly {
        path { meta, files -> "assembly/${meta.id}" }
    }
    assembly_logs {
        path { meta, files -> "assembly/${meta.id}/logs" }
    }
}

workflow {
    main:
    // ... workflow logic ...
    
    publish:
    assembly = ASSEMBLY_MAGS.out.assembly.map { meta, fasta -> [meta, fasta] }
                .mix(ASSEMBLY_MAGS.out.assembly_gfa)
    assembly_logs = ASSEMBLY_MAGS.out.assembly_log
}
```

---

## Conclusion

The Somatem pipeline is well-architected and follows many Nextflow best practices. The primary improvement opportunity is **migrating to modern global publishing syntax**, which will:

- Centralize output management
- Improve maintainability
- Enhance clarity of pipeline outputs
- Align with Nextflow v25.10+ best practices

Additional improvements in error handling, documentation, and testing will further strengthen the pipeline's robustness and usability.

**Recommended Next Steps:**
1. Implement global publishing refactor (see next document: REFACTOR_GUIDE.md)
2. Address high-priority items from section 11
3. Run `nextflow lint` and address all warnings
4. Test refactored pipeline with example datasets

---

**End of Code Review Report**
