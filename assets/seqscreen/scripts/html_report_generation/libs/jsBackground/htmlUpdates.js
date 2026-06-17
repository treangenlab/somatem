/* ######################################################################### */
/** Class to manipulate html */
var Html = function(columnCount) {
  /* #############################  INSTANCE VARIABLES  ###################### */
  this.columnCount = columnCount;
  /** Map sort icons to sort direction */
  this.sortIcon = {
    '1': 'fas fa-sort-up',
    '0': 'fas fa-sort',
    '-1': 'fas fa-sort-down'
  };

  /* #################################  FUNCTIONS  ########################### */

  /* ################################  TABLE STUFF  ########################## */

  /**
   * Draw table cells for all given sequences
   * @param {object} sequences sequences to be drawn
   */
  this.drawTable = function(sequences) {
    var res = '';
    var self = this;
    sequences.forEach(function(seq) {
      res +=
        '<tr>' +
        self.drawQueryCell(seq) +
        self.drawLengthCell(seq) +
        self.drawOrganismCell(seq) +
        self.drawGeneCell(seq) +
        self.drawGOCell(seq) +
        self.drawBsatCell(seq);
        for( var val in funsocs)
        {  
           res += self.drawFunsocsCell(seq[funsocs[val].name]) ;
        }
        res +='</tr>' +
        self.drawCollapsible(seq)+
        self.drawGOTermCollapsible(seq);
    });

    $('#results-table').html(res);
  };

  /** Draw query cell for sequence with optional link to alignment modal */
  this.drawQueryCell = function(seq) {
    var a = seq.query;
    if (seq.hasAlignment == 'true') {
      a =
        '<a href="#" data-toggle="modal" data-target="#alignment" data-seq="' +
        a + 
        '">' +
        a +
        '</a>';
    }

    return (
      '<td class="align-middle wrap" style="overflow-wrap: break-word;">' +
      a +
      '<span title="Show Sequence" class="faArrowIcon collapsed" ' +
      'data-toggle="collapse" data-target="#accordion-' +
      seq.qID +
      '"></span>'
    );
  };

  /** Draw length cell for sequence */
  this.drawLengthCell = function(seq) {
    return '<td class="align-middle centeredText wrap">' + seq.length + '</td>';
  };

  /** Draw organism cell for sequence */
  this.drawOrganismCell = function(seq) {
    return '<td class="align-middle centeredText wrap">' + seq.taxID + '</td>';
  };

  /** Draw gene cell for sequence */
  this.drawGeneCell = function(seq) {
    return '<td class="align-middle wrap">' + seq.geneName + '</td>';
  };

  /** Draw goTerm cell with first 3 go terms and link to visualization page*/
  this.drawGOCell = function(seq) {

    result =''
    if (seq.go_terms) {
      var goTerms = '<ol>';
      goTerms += seq.go_terms.slice(0, 3).reduce(function(res, term, idx) {
        return (res += '<li>'+ term.name + '</li>');
      }, '');
      goTerms +="</ol>"
      if(!(seq.go_terms.length==1 && ((seq.go_terms[0].id=="-")||(seq.go_terms[0].id==""))))
        result = '<a target="_blank"style="display:inline-block; width:95%" href="go.html?query=' +
        seq.query + '">' + goTerms +
        '</a>' + '<button onclick="changeText(this.id)" id="showMore'+seq.qID+'"type="button" title="Show GO-Term Breakdown" class="collapsed btn btn-light btn-xs" ' + 
        'data-toggle="collapse" data-target="#accordion-go-' +seq.qID+ 
        '">show more</button>'
      else 
        result = seq.go_terms[0].name
    }
    else result ='-'    
    
    return (
      '<td class="align-middle wrap" style="position:relative;">' +
      result+'</td>'
    );
  };
  /** Draw length cell for sequence */
  this.drawBsatCell = function (seq) {
    return '<td class="align-middle centeredText wrap">' + seq.bsat + '</td>';
  };

  /** Draw funsocs cell for sequence */
  this.drawFunsocsCell = function(funsocs) {
    var color;
    funsocs == '1' ? (color = '#ff99a1') : (color = null);
    return (
      '<td class="align-middle" >' +
      '<div style="min-width:20px; background:' +
      color +
      '; text-align:center; border-radius:20%">' +
      funsocs +
      '</div> </td>'
    );
  };

  /** Draw collapsible sequence cell for sequence */
  this.drawCollapsible = function(seq) {
    return (
      '<tr><td colspan="' +
      columnCount +
      '" class="hiddenRow bg-black ">' +
      '<div id="accordion-' +
      seq.qID +
      '" class="p-3 collapse collapseQuery">' +
      '</div></td></tr>'
    );
  };

  /** Fill collapsible sequence cell dynamically - is more efficient than filling all at once */
  this.fillCollapsible = function(seq) {
    var sequence = seq.sequence.split('').reduce(function(res, char) {
      return (res += '<span class="' + char + '">' + char + '</span>');
    }, '');

    var html =
      '<p class="text-break word-break text-monospace"' + sequence + '</p>';

    $('#accordion-' + seq.qID).html(html);
  };

  /** Draw collapsible confidence cell for sequence */
  this.drawGOTermCollapsible = function(seq) {
    var goTerms = '';
    if (seq.go_childs != null)
      seq.go_childs.forEach(function (goTerm) {
        goTerms +=
          '<span style="font-weight:bold">' +
          goTerm.id +
          '</span> - ' +
          goTerm.name +
          '<br>';
      });
    return '<tr><td colspan="' + columnCount + '" class="hiddenRow">' +
            '<div id="accordion-go-' + seq.qID + '" class="p-3 collapse"><h4> GO Terms</h4>' +                        
            goTerms + "</div></td></tr>";
    };

  /** Close all collapsibles */
  this.closeSecretRows = function() {
    $('collapse').collapse('hide');
  };

  /* ###############################  SORTING STUFF  ######################### */

  /** Reset all sortIcons to neutral state */
  this.resetSortIcons = function() {
    $('.sortTag').html('<i class="' + this.sortIcon[0] + '"></i>');
  };

  /** Set sortIcon at @param col to @param dir */
  this.setSortIcon = function(col, dir) {
    $('#columnResponse_'+col)[0].innerHTML = '<i class="' + this.sortIcon[dir] + '"></i>';
  };
  /* ######################################################################### */
  /** declaring variables */
     var checkBoxNames = [];
     var funsocs_data = funsocs;
  /* ######################################################################### */

  /* ################################  FILTER STUFF  ######################### */
  /**
  ** Assigning the FunSOCS names to CheckBoxnames
  **/
  this.assign_names=function(){
       var val = 0;
       var index = 5;
       for(val in funsocs_data)
       {
         checkBoxNames[index]= funsocs_data[val].name;
	 index++;
       }
   }
  /**
   * Draw all enabled filters as <ul> in #filterList div
   * @param {object} filters {col: content}
   */
  this.drawFilters = function(filters) {
    var res = '';
    var self = this;
    var clearButton ='';

    if (!$.isEmptyObject(filters)) {
      /** general open span element */
      var span = '<span class="badge">';
      // declare reduce functions here for readability
      res +=
        '<h4>Enabled Filters:' +
        '</h4><ul>';
      for (var f in filters) {
        res += '<li>';

        if (f == '0') {
          res += 'Query: ' + span + filters[f] + '</span>';
        }
        if (f == '1') {
          res += 'Sequence Length: ';
          if (filters[f].lower)
            res += span + 'Min ' + filters[f].lower + '</span> ';
          if (filters[f].upper)
            res += span + 'Max ' + filters[f].upper + '</span>';
        }
        if (f == '2') {
          res += 'Organism: ' + span + filters[f] + '</span>';
        }
        if (f == '3') {
          res += 'Gene Name: ' + span + filters[f] + '</span>';
        }
        if (f == '4') {
          res += 'GO Term: ' + span + filters[f] + '</span>';
        }
        if (parseInt(f) > 4 && parseInt(f) < columnCount) {
          res += checkBoxNames[parseInt(f)] + ': ' + span + ' </span>';
        }
        res += '</li>';
      }
      res += '</ul>';
      clearButton +=
        "<button id='clearAll' class='btn'>" +
        'Clear Filters ('+
        " <i class='fas fa-times'/>" +
        ' )'+
        ' </button>';
    }

    $('#filterList').html(res);
    $('#clearFilter').html(clearButton);
  };

  /**
   * Handle filter modal:
   *
   * Set inputs to old or empty value
   * Call reset() on clearButton
   * Get inputs and call success with non-empty inputs on applyButton
   */
  this.showFilter = function(prev, success, reset) {
    // Buttons
    var apply = $('#modalApply');
    var clear = $('#modalClear');
    // Filter contents
    var query = $('#modalQuery');
    var bounds = {
      lower: $('#modalMin'),
      upper: $('#modalMax')
    };
    var organism = $('#modalOrganism');
    var geneName = $('#modalGeneName');
    var goTerm = $('#modalGOTerm');

    var index = 0;
    var funSOCS = [];
    var minux = 5;
    for ( index =0; index < columnCount-minux; index++)
    {
      var name = "modalFunsocs" +index;
      funSOCS[index+minux] = document.getElementById(name);
    }
    // disable clearButton if no previous filter was set
    clear.prop('disabled', $.isEmptyObject(prev));

    // set current (old) filter values or empty
    query.val(prev[0] || '');
    if (!!prev[1]) {
      bounds.lower.val(prev[1].lower || '');
      bounds.upper.val(prev[1].upper || '');
    } else {
      bounds.lower.val('');
      bounds.upper.val('');
    }
    organism.val(prev[2] || '');
    geneName.val(prev[3] || '');
    goTerm.val(prev[4] || '');

    for (let i = 5; i < columnCount; i++) {
      var name = checkBoxNames[i];
      !!prev[i]
        ? funSOCS[i].checked =prev[i]
        : funSOCS[i].checked=false;
    }
    // clearButton functionality
    clear.off('click');
    clear.click(function() {
      reset();
    });

    // applyButton: send {col: content} as callback
    apply.off('click');
    apply.click(function() {
      // get all values, even if empty
      var contents = [ ] 
      var contents = {
        0: query.val(),
        1: { lower: bounds.lower.val(), upper: bounds.upper.val() },
        2: organism.val(),
        3: geneName.val(),
      };
      for(index = 5; index < columnCount; index++)
      {
        contents[index] = funSOCS[index].checked;
      }

      // remove falsy values (e.g. empty string, all checkboxes checked)
      for (let i = 0; i < columnCount; i++) {
        if (i == 1) {
          if (!contents[1].lower && !contents[1].upper) delete contents[1];
        } else if (!contents[i]) delete contents[i];
      }

      success(contents);
    });
  };

  /* ##############################  ALIGNMENT STUFF  ######################## */

  /** Fill alignment modal with sequence data
   * @param {object} seq selected Sequence
   * @param {$object} modal modal to apply this to
   */
  this.showAlignment = function(seq, modal) {
    $(modal)
      .find('.modal-title')
      .text('Alignment Results for ' + seq.query);

    $(modal)
      .find('#alignGene')
      .text(seq.geneName);
    $(modal)
      .find('#alignUpID')
      .text(seq.uniprot);
    $(modal)
      .find('#alignTaxo')
      .text(seq.taxID);
    //polyfill for "repeat"	  
    if (!String.prototype.repeat) {
        String.prototype.repeat = function(count) {
          'use strict';
          if (this == null)
            throw new TypeError('can\'t convert ' + this + ' to object');
      
          var str = '' + this;
          // To convert string to integer.
          count = +count;
          // Check NaN
          if (count != count)
            count = 0;
      
          if (count < 0)
            throw new RangeError('repeat count must be non-negative');
      
          if (count == Infinity)
            throw new RangeError('repeat count must be less than infinity');
      
          count = Math.floor(count);
          if (str.length == 0 || count == 0)
            return '';
      
          // Ensuring count is a 31-bit integer allows us to heavily optimize the
          // main part. But anyway, most current (August 2014) browsers can't handle
          // strings 1 << 28 chars or longer, so:
          if (str.length * count >= 1 << 28)
            throw new RangeError('repeat count must not overflow maximum string size');
      
          var maxCount = str.length * count;
          count = Math.floor(Math.log(count) / Math.log(2));
          while (count) {
             str += str;
             count--;
          }
          str += str.substring(0, maxCount - str.length);
          return str;
        }
    }
    var hsp = seq.alignment;
    var alignments = '';

    var qStart = parseInt(hsp.queryFrom);
    var hStart = parseInt(hsp.hitFrom);

    for (var index = 0; index < hsp.qSeq.length; index += 50) {
      var end = Math.min(hsp.qSeq.length, index + 50);

      qLine = hsp.qSeq.substr(index, 50);
      mLine = hsp.midLine.substr(index, 50).replace(/ /g, '&nbsp;');
      hLine = hsp.hSeq.substr(index, 50);

      qEnd = qStart + end - index - (qLine.split('-').length - 1); // count occurrences
      hEnd = hStart + end - index - (hLine.split('-').length - 1); // of '-' in string

      // sorry, this is as good as it gets
      alignments +=
        '<li class="list-group-item" style="font-family: monospace">' +
        '<p class="m-0 p-0">' +
        'Query' +
        '&nbsp;'.repeat(3) +
        qStart +
        '&nbsp;'.repeat(5 - qStart.toString().length) +
        qLine +
        '&nbsp;' +
        qEnd +
        '</p><p class="m-0 p-0">' +
        '&nbsp;'.repeat(13) +
        mLine +
        '</p><p class="m-0 p-0">' +
        'Subject&nbsp;' +
        hStart +
        '&nbsp;'.repeat(5 - hStart.toString().length) +
        hLine +
        '&nbsp;' +
        hEnd +
        '</p>' +
        '</li>';

      qStart = qEnd;
      hStart = hEnd;
    }

    $(modal)
      .find('#alignList')
      .html(alignments);
  };
};
