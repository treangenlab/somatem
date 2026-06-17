"use strict";

function _createForOfIteratorHelperLoose(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; return function () { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } it = o[Symbol.iterator](); return it.next.bind(it); }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

//Polyfill for Map function 
if (!Array.prototype.map) {

  Array.prototype.map = function(callback/*, thisArg*/) {

    var T, A, k;

    if (this == null) {
      throw new TypeError('this is null or not defined');
    }

    // 1. Let O be the result of calling ToObject passing the |this| 
    //    value as the argument.
    var O = Object(this);

    // 2. Let lenValue be the result of calling the Get internal 
    //    method of O with the argument "length".
    // 3. Let len be ToUint32(lenValue).
    var len = O.length >>> 0;

    // 4. If IsCallable(callback) is false, throw a TypeError exception.
    // See: http://es5.github.com/#x9.11
    if (typeof callback !== 'function') {
      throw new TypeError(callback + ' is not a function');
    }

    // 5. If thisArg was supplied, let T be thisArg; else let T be undefined.
    if (arguments.length > 1) {
      T = arguments[1];
    }

    // 6. Let A be a new array created as if by the expression new Array(len) 
    //    where Array is the standard built-in constructor with that name and 
    //    len is the value of len.
    A = new Array(len);

    // 7. Let k be 0
    k = 0;

    // 8. Repeat, while k < len
    while (k < len) {

      var kValue, mappedValue;

      // a. Let Pk be ToString(k).
      //   This is implicit for LHS operands of the in operator
      // b. Let kPresent be the result of calling the HasProperty internal 
      //    method of O with argument Pk.
      //   This step can be combined with c
      // c. If kPresent is true, then
      if (k in O) {

        // i. Let kValue be the result of calling the Get internal 
        //    method of O with argument Pk.
        kValue = O[k];

        // ii. Let mappedValue be the result of calling the Call internal 
        //     method of callback with T as the this value and argument 
        //     list containing kValue, k, and O.
        mappedValue = callback.call(T, kValue, k, O);

        // iii. Call the DefineOwnProperty internal method of A with arguments
        // Pk, Property Descriptor
        // { Value: mappedValue,
        //   Writable: true,
        //   Enumerable: true,
        //   Configurable: true },
        // and false.

        // In browsers that support Object.defineProperty, use the following:
        // Object.defineProperty(A, k, {
        //   value: mappedValue,
        //   writable: true,
        //   enumerable: true,
        //   configurable: true
        // });

        // For best browser support, use the following:
        A[k] = mappedValue;
      }
      // d. Increase k by 1.
      k++;
    }

    // 9. return A
    return A;
  };
}
//Polyfill for "find" function 
if (!Array.prototype.find) {
  Object.defineProperty(Array.prototype, 'find', {
    value: function(predicate) {
      // 1. Let O be ? ToObject(this value).
      if (this == null) {
        throw TypeError('"this" is null or not defined');
      }

      var o = Object(this);

      // 2. Let len be ? ToLength(? Get(O, "length")).
      var len = o.length >>> 0;

      // 3. If IsCallable(predicate) is false, throw a TypeError exception.
      if (typeof predicate !== 'function') {
        throw TypeError('predicate must be a function');
      }

      // 4. If thisArg was supplied, let T be thisArg; else let T be undefined.
      var thisArg = arguments[1];

      // 5. Let k be 0.
      var k = 0;

      // 6. Repeat, while k < len
      while (k < len) {
        // a. Let Pk be ! ToString(k).
        // b. Let kValue be ? Get(O, Pk).
        // c. Let testResult be ToBoolean(? Call(predicate, T, « kValue, k, O »)).
        // d. If testResult is true, return kValue.
        var kValue = o[k];
        if (predicate.call(thisArg, kValue, k, o)) {
          return kValue;
        }
        // e. Increase k by 1.
        k++;
      }

      // 7. Return undefined.
      return undefined;
    },
    configurable: true,
    writable: true
  });
}
//polyfill for "includes" function
if (!Array.prototype.includes) {
  Object.defineProperty(Array.prototype, 'includes', {
    value: function (searchElement, fromIndex) {
 
      if (this == null) {
        throw new TypeError('"this" is null or not defined');
      }
 
      // 1. Let O be ? ToObject(this value).
      var o = Object(this);
 
      // 2. Let len be ? ToLength(? Get(O, "length")).
      var len = o.length >>> 0;
 
      // 3. If len is 0, return false.
      if (len === 0) {
        return false;
      }
 
      // 4. Let n be ? ToInteger(fromIndex).
      //    (If fromIndex is undefined, this step produces the value 0.)
      var n = fromIndex | 0;
 
      // 5. If n ≥ 0, then
      //  a. Let k be n.
      // 6. Else n < 0,
      //  a. Let k be len + n.
      //  b. If k < 0, let k be 0.
      var k = Math.max(n >= 0 ? n : len - Math.abs(n), 0);
 
      function sameValueZero(x, y) {
        return x === y || (typeof x === 'number' && typeof y === 'number' && isNaN(x) && isNaN(y));
      }
 
      // 7. Repeat, while k < len
      while (k < len) {
        // a. Let elementK be the result of ? Get(O, ! ToString(k)).
        // b. If SameValueZero(searchElement, elementK) is true, return true.
        if (sameValueZero(o[k], searchElement)) {
          return true;
        }
        // c. Increase k by 1. 
        k++;
      }
 
      // 8. Return false
      return false;
    }
  });
}

/*Function to display loading symbol*/
document.onreadystatechange = function() {
  if (document.readyState !== "complete") {
      document.querySelector("#loader").style.visibility = "visible";
  } else {

      document.querySelector("#loader").style.visibility = "hidden";
      document.querySelector("body").style.visibility = "visible";
  }
};

//WIP
/// GLOBALS

/** data object used with the network */
var data = {
  nodes: new vis.DataSet(),
  edges: new vis.DataSet()
};
/** current sequence, found in all sequences using GET param 'query' */

var seq = sequences.find(function (seq) {
  return seq.query == new URLSearchParams(window.location.search).get('query');
});
/** ids of all go terms in seq */

var goIds = seq.go_terms.map(function (term) {
  return term.id;
});
/** ids of go terms without parents in seq */

var childsId = seq.go_childs.map(function (child) {
  return child.id;
});
/** relevant edges for this network */

var filteredNetwork = go_network.filter(function (item) {
  return goIds.includes(item.from) && goIds.includes(item.to);
});
/** all currently colored node ids, used by colorParents and resetColor*/

var coloredNodes = [];
/** set title */

document.getElementById('sequence').innerHTML = seq.query; /// helper functions

/**
 * Wrap text after max. limit characters, don't break words
 * @param {string} text 
 * @param {number} limit 
 */

var wordwrap = function wordwrap(text, limit) {
  var lines = [];
  var current = '';

  for (var _iterator = _createForOfIteratorHelperLoose(text.split(' ')), _step; !(_step = _iterator()).done;) {
    var word = _step.value;

    if (current.length + word.length > limit) {
      lines.push(current);
      current = word;
    } else {
      current += " " + word;
    }
  }

  lines.push(current); // push last line 

  return lines.join("\n").trim(); // first line starts with space
};
/**
 * return rgb color code for node at index
 * @param {object} node 
 * @param {number} index # of node's go term in seq.go_terms
 */


var getColor = function getColor(node, index) {
  if (['molecular_function', 'biological_process', 'cellular_component'].includes(node.label || node.name)) return '#007bff';
  if (childsId.includes(node.id)) return '#28a745';
  return '#17a2b8'; // ???, bootstrap info
};

var getLevel = function getLevel(term, filter, result) {
  return _getLevel(term.id, filter, result);
};

var _getLevel = function _getLevel(id, filter, result) {
  // return 0 for top-level terms
  if (['GO:0003674', 'GO:0008150', 'GO:0005575'].includes(id)) return 0; // find parents, return level + 1 if parents already in results → much faster

  var parents = filter.filter(function (item) {
    return item.to == id;
  }).map(function (parent) {
    return parent.from;
  });
  console.log(result, parents, id);
  if (parents.every(function (parent) {
    return parent in result;
  })) return Math.max.apply(Math, parents.map(function (parent) {
    return result[parent];
  })) + 1;
  return Math.max.apply(Math, parents.map(function (parent) {
    return _getLevel(parent, filter, result) + 1;
  }));
};
/**
 * TBD: try to get these things done in a readable and efficient way:
 * * calculate each node exactly once
 * * if calculating parent nodes is necessary, save nodes and parent nodes (or return? {${id}: [level, parents]})
 * @param {number} id id of node to be calculated
 * @param {object} result {id: {level, parents}}
 * @return {object} {id: {level, parents}}
 */


var getStructData = function getStructData(id, filter, result) {
  // return result if id has been calculated before
  if (id in result) return result; // return 0, no parents for top-level terms

  if (['GO:0003674', 'GO:0008150', 'GO:0005575'].includes(id)) {
    result[id] = {
      level: 0,
      parents: []
    };
    return result;
  } // find parents, return level + 1 if parents already in results → much faster


  var parents = filter.filter(function (item) {
    return item.to == id;
  }).map(function (parent) {
    return parent.from;
  });

  if (parents.every(function (parent) {
    return parent in result;
  })) {
    var _level = Math.max.apply(Math, parents.map(function (parent) {
      return result[parent].level;
    })) + 1;

    if (_level < 0) _level = 0;
    result[id] = {
      level: _level,
      parents: parents
    };
    return result;
  }

  var level = Math.max.apply(Math, parents.map(function (parent) {
    return getStructData(parent, filter, result)[parent].level + 1;
  }));
  if (level < 0) level = 0;
  result[id] = {
    level: level,
    parents: parents
  };
  return result;
};
/**
 * color a node and all its parent nodes
 */


var colorParents = function colorParents(id) {
  // can safely assume that a node will be found
  var node = data.nodes.get().find(function (node) {
    return node.id === id;
  });
  /* if (!node) {
      node = go.find(node => node.id == id)
       data.nodes.add(node)
  } */

  coloredNodes.push(node);
  var parents = node.parents;

  if (parents.length) {
    parents.forEach(function (parent) {
      return colorParents(parent);
    });
  } else {
    coloredNodes.forEach(function (node) {
      return node.color = '#e5da2a';
    });
    data.nodes.update(coloredNodes);
  }
};
/**reset colorParents(id) */


var resetColor = function resetColor(_) {
  coloredNodes.forEach(function (node) {
    return node.color = getColor(node, node.index);
  });
  data.nodes.update(coloredNodes);
  coloredNodes = [];
}; // calculate hierarchy level for each node/go term


var meta = {};
var level = {}; // map go_terms to nodes with calculated level and coordinates

var go = seq.go_terms.map(function (term, index) {
  if (!(term.id in meta)) {
    // don't do anything if this term has been calculated before
    var result = getStructData(term.id, filteredNetwork, meta); // get data of term and all its parents

    for (var id in result) {
      meta[id] = result[id]; // add id to level object

      var l = result[id].level;
      if (!level[l]) level[l] = []; // if id's level doesn't yet exist in level-object add empty array

      level[l].push(id);
    }
  }

  return {
    id: term.id,
    index: index,
    label: wordwrap(term.name, 20),
    title: term.id + " " + term.name,
    level: meta[term.id].level,
    parents: meta[term.id].parents,
    value: 2,
    color: getColor(term, index)
  };
});
var nodes = [];

for (var l in level) {
  nodes.push.apply(nodes, new Set(level[l]));
} //nodes.forEach((node, index) => setTimeout(_ => data.nodes.add(go.find(n => n.id == node)), index * 50))
// set data
//data.nodes = new vis.DataSet(go.filter(node => !filteredNetwork.find(edge => edge.from == node.id)))


data.nodes = new vis.DataSet(go);
data.edges = new vis.DataSet(filteredNetwork);
var options = {
  nodes: {
    shape: 'dot',
    scaling: {
      min: 20,
      max: 30,
      label: {
        min: 14,
        max: 30,
        drawThreshold: 9,
        maxVisible: 20
      }
    },
    font: {
      size: 14,
      face: 'Helvetica Neue, Helvetica, Arial'
    }
  },
  layout: {
    hierarchical: {
      direction: 'UD',
      levelSeparation: 250,
      nodeSpacing: 150,
      treeSpacing: 500,
      sortMethod: 'directed'
    }
  },
  interaction: {
    hover: true,
    hoverConnectedEdges: false,
    selectConnectedEdges: true
  },
  physics: {
    hierarchicalRepulsion: {
      nodeDistance: 250
    },
    stabilization: true
  }
};
var container = document.getElementById('network');
var network = new vis.Network(container, data, options);
/*
network.on('click', args => {
if (!args.nodes.length) return
let node = go.find(node => node.id == args.nodes[0])
let parents = go.filter(n => node.parents.includes(n.id))
let allNodes = data.nodes.get()
allNodes.forEach(node => node.hidden = true)
parents.forEach(parent => {
parent.hidden = false
})
data.nodes.update(allNodes)
data.nodes.update(parents)
})*/

network.on('hoverNode', function (args) {
  return colorParents(args.node);
});
network.on('blurNode', function (_) {
  return resetColor();
});
