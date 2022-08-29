function show_hide_params(){
  var l = document.getElementById("rlinkage");
  var d = document.getElementById("rdistance");
  var n = document.getElementById("rnclusters");
  var e = document.getElementById("algorithm");  
  l.style.display = "none";
  d.style.display = "none";
  n.style.display = "none";
    if(e.value=='hierarchical'){
	      l.style.display = "table-row";
    }
    else if(e.value=='kmeans'){
          n.style.display = "table-row";
    }
    else if(e.value=='tda'){
      l.style.display = "none";
      d.style.display = "none";
      n.style.display = "none";
    }
    else{
      d.style.display = "table-row";
    }   
}

function download_text(elementId,ext='text') {
    var text = document.getElementById(elementId).innerHTML;
    text=text.replaceAll("&amp;","&");
    var link = document.createElement('a');
    alert(text);
    link.setAttribute('download', 'table.'+ext);
    link.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    link.click(); 
}

function show_hide_params_rd(){
  var o = document.getElementById("tr_order");
  var t = document.getElementById("tr_threshold");
  var e = document.getElementById("algorithm");  
  o.style.display = "none";
  t.style.display = "none";  
  if(e.value == 0 ){
	      o.style.display = "table-row";
  }
  else{
      t.style.display = "table-row";
  }   
}
