/* ######################################################################### */
(function() {
  'use strict';

  /* ##################################  GLOBALS  ############################ */

  /** current sorting state */
  var sort = {
    column: 0,
    direction: 1 // 1 or -1
  };
  /** Variable for pagination */ 
  var currentPage = 1;
  var numberPerPage = 100;
  var numberOfPages = 1;
  var numberOfButtons = 8;
  var fileread_flag = 0;
  var numFiles = 0;

  /** functions to sort each column by */
  // Regular expression to separate the digit string from the non-digit strings.
  var reParts = /\d+|\D+/g;
 
  // Regular expression to test if the string has a digit.
  var reDigit = /\d/;
 
  // Add cmpStringsWithNumbers to the global namespace.  This function takes to
  // strings and compares them, returning -1 if `a` comes before `b`, 0 if `a`
  // and `b` are equal, and 1 if `a` comes after `b`.
  function cmpStringsWithNumbers(a, b) {
    // Get rid of casing issues.
    a = a.toUpperCase();
    b = b.toUpperCase();
    
    //"-" sorts after everything else
    if (a =="-") return -1;
    if (b =="-") return 1;
  
    // Separates the strings into substrings that have only digits and those
    // that have no digits.
    var aParts = a.match(reParts);
    var bParts = b.match(reParts);
 
    // Used to determine if aPart and bPart are digits.
    var isDigitPart;
 
    // If `a` and `b` are strings with substring parts that match...
    if(aParts && bParts &&
        (isDigitPart = reDigit.test(aParts[0])) == reDigit.test(bParts[0])) {
      // Loop through each substring part to compare the overall strings.
      var len = Math.min(aParts.length, bParts.length);
      for(var i = 0; i < len; i++) {
        var aPart = aParts[i];
        var bPart = bParts[i];
 
        // If comparing digits, convert them to numbers (assuming base 10).
        if(isDigitPart) {
          aPart = parseInt(aPart, 10);
          bPart = parseInt(bPart, 10);
        }
 
        // If the substrings aren't equal, return either -1 or 1.
        if(aPart != bPart) {
          return aPart < bPart ? 1 : -1;
        }
 
        // Toggle the value of isDigitPart since the parts will alternate.
        isDigitPart = !isDigitPart;
      }
    }
 
    // Use normal comparison.
    return (a >= b) - (a <= b);
  };
  // This function takes FUNCOS as argument and compares them,
  // Expected order is "1,0,NA,-"
  function cmpfuncos(a,b){
    //"-" sorts after everything else
    if (a =="-") return -1;
    if (b =="-") return 1;
    
    //"NA" sorts after "-"
    if (a =="NA") return -1;
    if (b =="NA") return 1;
 
    return a > b ? 1 : a < b ? -1 : 0;
   }
   function sortFuncos(index, a, b)
   {
     var col_name = funsocs[index].name;
     return cmpfuncos(a[col_name],b[col_name]);
   }
  var sortFunctions = {
    0: function(a, b) {
      return cmpStringsWithNumbers(a.query, b.query);
    },
    1: function(a, b) {
      return a.length - b.length;
    },
    2: function(a, b) {
      return cmpStringsWithNumbers(a.taxID, b.taxID);
    },
    3: function(a, b) {
      return cmpStringsWithNumbers(a.geneName, b.geneName);
    }
  };
  /** active filters */
  var filter = {};
  /** functions to filter each column by */
  var filterFunctions = {
    0: function(item) {
      return item.toLowerCase().includes(filter[0].toLowerCase());
    },
    1: function(item) {
      return (
        (!filter[1].lower || parseInt(filter[1].lower) <= item) &&
        (!filter[1].upper || parseInt(filter[1].upper) >= item)
      );
    },
    2: function(item) {
      return item.toLowerCase().includes(filter[2].toLowerCase());
    },
    3: function(item) {
      return item.toLowerCase().includes(filter[3].toLowerCase());
    },
    4: function(item) {
      return item.toLowerCase().includes(filter[4].toLowerCase());
    },
    5: function(item) {
      return parseInt(item) == true;
    }
  };
  /** map column number to name */
  var columnNames = {
    0: 'query',
    1: 'length',
    2: 'taxID',
    3: 'geneName',
    4: 'goTerm'
  };
  
  var index = 5;
  for( var val in funsocs)
  {
       columnNames[index]= funsocs[val].name;
       index++;
  }
  /** number of columns */
  var threat_counts = [0, 0, 0, 0, 0];
  var columnCount = column_count;
  var arraySequences =[];
  var currentSequences = {};
  var validSequences = [];
  var sequences={};
  var sequence_length = 0; 
  /** functions to manipulate html */
  var html = new Html(columnCount);

  /* #############################  INITIAL SETUP  ########################### */

  $(document).ready(function() {
    // Preset sorting
    html.resetSortIcons();
    html.assign_names();
   
    // Upon File upload this function is called
    let fileUpload = document.getElementById("inputfile");
    fileUpload.onchange = function(event) {
      numFiles = event.target.files.length;

      //Restoring default values so that only new values are uploaded
      fileread_flag = 0;
      arraySequences =[];
      currentSequences = {};
      validSequences = [];
      sequences={};
      sequence_length = 0;
      document.getElementById('upload_msg').innerHTML = "File reading on progress";
      var i;
      for (i = 0; i < numFiles; i++) {
        var name = event.target.files[i].name;
        if(name.match(/output_part/g)== null)
        {
            document.getElementById('upload_msg').innerHTML = "unexpected file uploaded";
            document.getElementById("upload_msg").classList.add('btn-outline-danger');
            window.alert("Unexpected file uploaded");
            loadSortedTable(0);
            display_default_summary();
            return;
        }
      }
      let sortedFiles = [];
      for (let i = 0; i < numFiles; i++) {
        sortedFiles.push([i, event.target.files[i].name]);
      }
      sortedFiles.sort(function (a, b) {
        let compare = 0;
        if (b[1] < a[1]) {
          compare = 1;
        } else if (a[1] < b[1]) {
          compare = -1;
        }
        return compare;
      });

      for (i = 0; i < numFiles; i++) {
        var reader = new FileReader();
        reader.onload = onReaderLoad;
        reader.readAsText(event.target.files[i]);
      }
    }


    // Add click handlers
    $('#clearAll').click(function() {
      resetFilters();
    });
    validateFilterInput();
    
    $('.sortTag').click( function() {
      var num = Number($(this).attr('value'));
      loadSortedTable(num);
    });

    // Add handlers to modals
    $('#modalButton').click(function() {
      showFilter();
    });
    // function to pass 'this' (alignment-modal)
    $('#alignment').on('show.bs.modal', function(event) {
      showAlignment(event, this);
    });
  });

  /* ################################  FILTER ERROR HANDLING  ########################### */
  function validateFilterInput(){
    var seqmin = document.getElementById('modalMin');
    var seqmax = document.getElementById('modalMax');
    seqmin.addEventListener('input',function(){
        var seqmin = this.value;
        var seqmax = document.getElementById('modalMax').value;
        var minval = Number(seqmin);
        var maxval = Number(seqmax);
        if((minval < 1) &(seqmin != ""))
        {
            document.getElementById('sequence_warning').innerHTML = "Mininum Value is 1";            
            document.getElementById("modalApply").disabled = true;
        }
        else if((maxval < 1)&(seqmax != ""))
        {
            document.getElementById('sequence_warning').innerHTML = "Mininum Value is 1";            
            document.getElementById("modalApply").disabled = true; 
        }
        else if((minval > maxval)&(seqmax != "")&(seqmin != ""))
        {
            document.getElementById('sequence_warning').innerHTML = "Mininum value must be less than the maximum value";            
            document.getElementById("modalApply").disabled = true; 
        }
        else
        {
            document.getElementById('sequence_warning').innerHTML = "";            
            document.getElementById("modalApply").disabled = false; 
        }
     }, false); 
    seqmax.addEventListener('input',function(){
        var seqmin = document.getElementById('modalMin').value;
        var seqmax = this.value;
        var minval = Number(seqmin);
        var maxval = Number(seqmax);
        if((minval < 1) &(seqmin != ""))
        {
            document.getElementById('sequence_warning').innerHTML = "Mininum Value is 1";            
            document.getElementById("modalApply").disabled = true;
        }
        else if((maxval < 1)&(seqmax != ""))
        {
            document.getElementById('sequence_warning').innerHTML = "Mininum Value is 1";            
            document.getElementById("modalApply").disabled = true; 
        }
        else if((minval > maxval)&(seqmax != "")&(seqmin != ""))
        {
            document.getElementById('sequence_warning').innerHTML = "Mininum value must be less than the maximum value";            
            document.getElementById("modalApply").disabled = true; 
        }
        else
        {
            document.getElementById('sequence_warning').innerHTML = "";            
            document.getElementById("modalApply").disabled = false; 
        }
     }, false);
  }

  //Function to read the files uploaded
  function onReaderLoad(event){

    var tmpSequences = JSON.parse(event.target.result);
    display_summary(tmpSequences.shift());

        var tmpArraySeqs = Object.keys(tmpSequences).map(function (key) {
            return tmpSequences[key];
        });

        var tmpOldArraySeqs = Object.keys(sequences).map(function (key) {
            return sequences[key];
        });

        var newSeqs = tmpOldArraySeqs.concat(tmpArraySeqs);
        sequences = { ...newSeqs };

        var tmpArraySequences = Object.keys(tmpSequences).map(function (key) {
            return tmpSequences[key];
        });
        fileread_flag += 1;
        arraySequences = arraySequences.concat(tmpArraySequences);
        validSequences = validSequences.concat(tmpArraySequences);

        if (fileread_flag == numFiles) {
            loadSortedTable(0);
            document.getElementById('upload_msg').innerHTML = "File reading completed";
        }
  }
    
  //Function that displays the summary information
  function display_summary(summary)
  {
    sequence_length = sequence_length + summary.length;
    document.getElementById('mode').innerHTML = summary.mode;
    document.getElementById('flag').innerHTML = "Runtime Flag: "+summary.rflag;
    document.getElementById('fasta_file').innerHTML = summary.fasta_file;
    document.getElementById('input_file').innerHTML = '';
    var file;
    for(file in summary.mtimes)
    {
      if(file == "tsv")
      {
        document.getElementById('input_file').innerHTML +=' <div class="badge badge-primary"> S2FAST report '+ summary.mtimes[file]+"</div>"; 
      }
      else if (file == "fasta")
      {
        document.getElementById('input_file').innerHTML +=  ' <div class="badge badge-primary"> Sequence file '+ summary.mtimes[file]+'</div>'; 
      }
      else if (file == "xml")
      {
        document.getElementById('input_file').innerHTML +=  ' <div class="badge badge-primary"> Blast output '+ summary.mtimes[file]+'</div>'; 
      }
    }
  document.getElementById('seq_len').innerHTML = "Sequences Provided: " + sequence_length;
  }

  //Function that displays the default summary information
function display_default_summary()
  {
    document.getElementById('mode').innerHTML = "";

    document.getElementById('flag').innerHTML = "";
    document.getElementById('fasta_file').innerHTML = "";
    document.getElementById('seq_len').innerHTML = "";
    document.getElementById('input_file').innerHTML =""; 
 }

 
  /* ################################  PAGINATION  ########################### */
  /**
   * Function to initialize paging
   */
  function initializePaging(){
    numberOfPages = Math.ceil(validSequences.length / numberPerPage);
    currentPage = 1;
    if(numberOfPages > 1)
    {
       var start = (currentPage -1)*numberPerPage;
       var end = currentPage * numberPerPage;
       currentSequences = validSequences.slice(start, end);
       pageButtons();
    }
    else
    {
      currentSequences = validSequences;
      $('#pagination-wrapper').empty();
      if(validSequences == 0)
      {
        var wrapper = document.getElementById('pagination-wrapper');
        wrapper.innerHTML = 'No Records Found';
        wrapper = document.getElementById('pagination-wrapper-top');
        wrapper.innerHTML = 'No Records Found';
      }
      else
      {
        var wrapper = document.getElementById('pagination-wrapper');
        var page = 1;
        var str = "Top"+page;
        wrapper.innerHTML = '<button value='+page+' id='+page+' class="page btn btn-sm btn-info">'+page+'</button>';
        wrapper = document.getElementById('pagination-wrapper-top');
        wrapper.innerHTML = '<button value='+page+' id='+str+' class="page btn btn-sm btn-info">'+page+'</button>';
	document.getElementById(page).disabled = true;
        document.getElementById(str).disabled = true;
      }
    }
  }
  /**
   * Function to update page buttons
   */
  function pageButtons() {
    var wrapper = document.getElementById('pagination-wrapper');
    var wrapperTop = document.getElementById('pagination-wrapper-top');
    wrapper.innerHTML = '';
    wrapperTop.innerHTML = '';
    var maxLeft = (currentPage - Math.floor(numberOfButtons / 2));
    var maxRight = (currentPage + Math.floor(numberOfButtons/ 2));

    if (maxLeft < 1) {
        maxLeft = 1;
        maxRight = numberOfButtons;
    }

    if (maxRight > numberOfPages) {
        maxLeft = numberOfPages - (numberOfButtons - 1);

        if (maxLeft < 1){
        	maxLeft = 1;
        }
        maxRight = numberOfPages;
    }
    for (var page = maxLeft; page <= maxRight; page++) {
        var str = "Top"+page;
	wrapper.innerHTML += '<button value='+page+' id='+page+' class="page btn btn-sm btn-info">'+page+'</button>';
        wrapperTop.innerHTML += '<button value='+page+' id='+str+' class="page btn btn-sm btn-info">'+page+'</button>';
    }

    if (currentPage != 1) {
        wrapper.innerHTML = '<button value=0 id="prev" class="page btn btn-sm btn-info"> Prev </button>' + wrapper.innerHTML;
        wrapperTop.innerHTML = '<button value=0 id="prevTop" class="page btn btn-sm btn-info">Prev </button>' + wrapperTop.innerHTML;

    }
    if(maxLeft >= 2){
        wrapper.innerHTML = '<button value=1 class="page btn btn-sm btn-info"> First</button>' + wrapper.innerHTML;
        wrapperTop.innerHTML = '<button value=1 class="page btn btn-sm btn-info">First</button>' + wrapperTop.innerHTML;
    }
    if (currentPage != numberOfPages) {
        wrapper.innerHTML += '<button value=0 id="next" class="page btn btn-sm btn-info">Next</button>';
        wrapperTop.innerHTML += '<button value=0 id="nextTop" class="page btn btn-sm btn-info">Next</button>';
    }

    if(maxRight != numberOfPages){
        wrapper.innerHTML += '<button value='+numberOfPages+' class="page btn btn-sm btn-info">Last</button>';
        wrapperTop.innerHTML += '<button value='+numberOfPages+' class="page btn btn-sm btn-info">Last</button>';
    }
    var str = "Top"+currentPage;
    document.getElementById(currentPage).disabled = true;
    document.getElementById(str).disabled = true;
    
    $('.page').click( function() {
        $('#results-table').empty();
        var num = Number($(this).val());
	if(num == 0)
	{ 
	   var name = (String(this.id)).trim();
		
           if(name.valueOf() == new String("prev").valueOf())
	   {
	      currentPage = currentPage -1;
	   }
	   else if(name.valueOf() == new String("next").valueOf())
	   {
	      currentPage = currentPage +1;
	   }
	   else if(name.valueOf() == new String("prevTop").valueOf())
	   {
	      currentPage = currentPage -1;
	   }
	   else if(name.valueOf() == new String("nextTop").valueOf())
	   {
	      currentPage = currentPage +1;
	   }
	
	}
	else
	{
           currentPage = num;
	}
	updatePage();
       document.body.scrollTop = 350; // For Safari
       document.documentElement.scrollTop = 350; // For Chrome, Firefox, IE and Opera
    });

  }
  /**
   * Function to update page when corresponding page button is clicked
   */
  function updatePage(){
    var start = (currentPage -1)*numberPerPage;
    var end = currentPage * numberPerPage;
    currentSequences = validSequences.slice(start, end);
    html.drawTable(currentSequences);
    html.resetSortIcons();
    html.setSortIcon(sort.column, sort.direction);
    addCollapseFill();
    pageButtons();
  }
  
  /* ################################  SHOW STUFF  ########################### */

  /**
   * Fill alignment modal with data
   * @param {*} event.relatedTarget whatever triggered the modal
   * @param {object} modal modal to be triggered
   */
  function showAlignment(event, modal) {
    var value = $(event.relatedTarget).data('seq');
    var seq = arraySequences.find(function (seq) {
        return seq.query==value;
    });
    html.showAlignment(seq, modal);
  }

  /**
   * Fill collapse with data
   *
   * Useful because coloring each character is expensive
   * @param {*} event.relatedTarget whatever triggered the modal
   */
  function showCollapse(event) {
    var value = parseInt($(event.currentTarget).attr('id').split('-')[1]);
    var seq = arraySequences.find(function (seq) {
      return seq.qID==value;
      });
    html.fillCollapsible(seq);
  }

  /**
   * Handle callbacks for 'All Filters' modal
   * No params
   */
  function showFilter() {
    function success(newFilters) {
      resetFilters();

      for (var index in newFilters) {
        addFilter(index, newFilters[index]);
      }
      applyFilters();
    }

    function reset() {
      resetFilters();
    }

    html.showFilter(filter, success, reset);
  }

  /* ############################  FILTER / SORT STUFF  ###################### */

  function loadSortedTable(col) {
    if (validSequences.length < 50) sortTable(col); //no loading screen
    else {
      document.getElementById('cover-spin').style.visibility = 'visible';
      setTimeout(function() {
        sortTable(col).then(function() {
          document.getElementById('cover-spin').style.visibility = 'hidden';
        });
      }, 100);
    }
  }
  /**
   * Sort table by column
   * If table is already sorted by given column, inverse order
   * @param {number|string} col
   */
  function sortTable(col) {
    return new Promise(function(resolve) {
      if (sort.column === col) {
           sort.direction *= -1;
	}
       else {
        sort = {
          column: col,
          direction: -1
        };
      }

      validSequences.sort(function(a, b) {
          if(sort.column > 4)
             return sortFuncos(sort.column - 5,a,b) *sort.direction;
          else
             return sortFunctions[sort.column](a, b) * sort.direction;
       
      });
      
      //Add Pagination
      initializePaging();
      html.drawTable(currentSequences);
      html.resetSortIcons();
      html.setSortIcon(sort.column, sort.direction);
      addCollapseFill();
      setTimeout(function() {
        resolve();
      }, 0);
    });
  }

  /**
   * Add filter at column to filter-object
   * @param {number|string} col
   * @param {string} newFilter
   */
  function addFilter(col, newFilter) {
    filter[col] = newFilter;
  }

  /**
   * Apply all filters to sequences and store in currentSequences
   * No params
   */
  function applyFilters() {
    validSequences = arraySequences.filter(function(seq) {
      var accept = true;

      for (var index in filter) {
        // apply filterFunction at index with (column at index) of sequence
        // e.g. apply filterFunction[1] with seq[columnNames[1]] (→ seq['length'])
        var filterindex = index;
        if (index > 5 && index < columnCount) filterindex = 5; //use same filterfunction for checkboxes
        accept =
          accept && filterFunctions[filterindex](seq[columnNames[index]]);
      }

      return accept;
    });
    //Add Pagination
    initializePaging();
    html.drawTable(currentSequences);
    html.drawFilters(filter);
    addCollapseFill();
    $(document).ready(function() {
      // Add click handlers
      $('#clearAll').click(function() {
        resetFilters();
      });
    });
  }

  /**
   * Remove all filters, reset to initial state
   * No params
   */
  function resetFilters() {
    filter = {};
    applyFilters();
  }

  /**
   * (re-) add showCollapse to show.bs.collapse event for all .collapseQuery.
   * To be used everytime the table is rebuilt
   */
  function addCollapseFill() {
    $('.collapseQuery').off('show.bs.collapse');
    $('.collapseQuery').on('show.bs.collapse', function(event) {
      showCollapse(event);
    });
  }
})(); // anonymous wrapper

/**
 * Change text of the show-more-button on click
 */
function changeText(id){
  var buttonText = document.getElementById(id).innerHTML
  if(buttonText=="show more")
    document.getElementById(id).innerHTML="show less";
  
  else
  document.getElementById(id).innerHTML = "show more"
}
