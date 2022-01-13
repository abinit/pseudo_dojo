/*
function.prototype.method = function (name func) {
    this.prototype[name] = func
    return this
}
*/

var keys = ['hh', 'hl', 'hn', 'nv', 'd', 'dp', 'gb'];

var els = [
  'H', 'He', 
  'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne','Na', 'Mg', "Al", "Si", 'P', 'S', 'Cl', 'Ar', 
  'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
  'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
  'Cs','Ba',
  'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
  'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'
];

function getParameterByName(name) {
    var url = window.location.href;
    console.log(url);
    var name = name.replace(/[\[\]]/g, "\\$&");
    var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)");
    var results = regex.exec(url);
    if (!results) return null;
    if (!results[2]) return '';
    return decodeURIComponent(results[2].replace(/\+/g, " "));
}

function set_warning(txt) {
  // Set the text in the warning box
  var warningbox = document.getElementById('warning_box');
  warningbox.innerHTML = "<div class='alert warning'><span id='cbn' class='closebtn'>&times;</span><strong>Warning!</strong> ".concat(txt, "</div>");
  var close = document.getElementById("cbn");
  close.onclick = function(){
     var div = document.getElementById('warning_box');
     setTimeout(function(){ div.innerHTML = ""; }, 100);
  }
}


function make_light() {
    document.getElementById('FMT').value = 'psp8'
    hide_clases = ["hide", "name", 'intro', "styled-longselect", "selection_bar", "help_button", "description", "menubar"];
    for (var j = 0; j < hide_clases.length; j++) {
        var tohide = document.getElementsByClassName(hide_clases[j]);
        for (var i =0; i < tohide.length; i++) {
            tohide[i].style.visibility="hidden";
        }
    }
    document.getElementById('X_n').setAttribute("style","left:326px; top:91px; height:170px; width:140px;");
    document.getElementById('N').setAttribute("style","left:326px; top:91px; height:170px; width:140px; font-size=20px");
    document.getElementById("download_button").setAttribute("style","left:70px; top:151px; width:200px; height:55px; padding:15px");
    elements = document.getElementsByClassName('element')
    for (var i; i < elements.length; i++){
       elements[i].setAttribute('style', 'font-size:24px; margin-top:12px; line-height:1; text-align:center;');
    }
    document.getElementById("X_el").setAttribute('style', 'margin-top:20px;');
    document.getElementById("X_hl").setAttribute('style', 'font-size:20px; padding:2px');
    document.getElementById("X_hn").setAttribute('style', 'font-size:20px; padding:2px');
    document.getElementById("X_hh").setAttribute('style', 'font-size:20px; padding:2px');
    document.getElementById("X_nv").setAttribute('style', 'font-size:20px; margin-top:-158px; padding:2px');
    document.getElementById("det_test").setAttribute('style', 'font-size:20px; padding:2px');
    document.getElementById("det_hints").setAttribute('style', 'font-size:20px; margin-top:5px; padding:2px');
    document.getElementById("X_d").setAttribute('style', 'font-size:20px; padding:2px');
    document.getElementById("X_dp").setAttribute('style', 'font-size:20px; padding:2px');
    document.getElementById("X_gb").setAttribute('style', 'font-size:20px; padding:2px');
}

function set_info(info, animate) {
    var averages = [0,0,0,0,0,0,0];
    var sums = [0,0,0,0,0,0,0];
    if (animate * localStorage.getItem('animate') === 1){
        console.log('added animating');
        $('.plugin').removeClass('anim');
        $('.plugin').removeClass('chaos');
        setTimeout("$('.plugin').addClass('anim')",10)
    }
    for (el in els) {
        for (key in keys){
            var id_key = els[el] + '_' + keys[key];
            var x = document.getElementById(id_key);
            var el_info = info[els[el]];
            var val = 'na';
            if (typeof(el_info) == 'undefined') {
                val = 'na';
            }
            else {
                val = el_info[keys[key]];
            }
            if (val === 'na' || val === 'nan'){
                var xx = 1;
            }
            else {
                averages[key] += parseFloat(val);
                sums[key] += 1;
            }
            x.innerHTML = val;
        }
    }
    for (i in averages){
        averages[i] = averages[i]/sums[i]
        averages[i] = averages[i].toFixed(2)
    }
    var summary = "The averages of this table:\n\n"
    summary += "low hint\t\t\t\t: " + averages[1] + " Ha\n"
    summary += "normal hint\t\t\t: " + averages[2] + " Ha\n"
    summary += "high hint\t\t\t\t: " + averages[0] + " Ha\n"
    summary += "nuber of valence shells\t: " + averages[3] + " \n"
    summary += "delta\t\t\t\t\t: " + averages[4] + " meV\n"
    summary += "normalized delta\t\t: " + averages[5] + " \n"
    summary += "gbrv\t\t\t\t\t: " + averages[6] + " %\n"
    if (sums[0] > 0) {
        set_av(averages);
        reset_X()
    }
}

function loadJSON(file, callback) {
    // Helper function to load a json file.
    var xobj = new XMLHttpRequest();
    xobj.overrideMimeType("application/json");
    xobj.open('GET', file, true);
    xobj.onreadystatechange = function () {
        if (xobj.readyState == 4 && xobj.status == 200) {
            // Required use of an anonymous callback as .open will NOT return a value but simply returns undefined in asynchronous mode
            callback(xobj.responseText);
          }
    };
    xobj.send(null);
 }

function store_available_files() {
    // Get list of files from files.json and store it in the localStorage.
    loadJSON('files.json', function(response) {
        var info = JSON.parse(response);
        localStorage.setItem('files', info);
    });
}

function load_set_info(animate) {
    acc = document.getElementById('ACC').value;
    xcf = document.getElementById('XCF').value;
    type = document.getElementById('TYP').value;
    set_info({}, 0);
    var file = type + '_' + xcf + '_' + acc + '.json';
    loadJSON(file, function(response) {
        var info = JSON.parse(response);
        set_info(info, animate);
    });
}

function set_X(elm, color, n){
    // Update params shown in the X_n box.
    if (els.indexOf(elm) >= 0){
        document.getElementById('N').innerHTML = n;
        var x = document.getElementById('X_n');
        x.style.backgroundColor = color;
        var x = document.getElementById('X_el');
        x.innerHTML = elm;
        for (key in keys){
            var id_key = 'X_' + keys[key];
            var id_key_in = elm + '_' + keys[key];
            // Get the params from the pseudo associated to this element and copy to the X box.
            var x = document.getElementById(id_key_in);
            var y = document.getElementById(id_key);
            var val = x.innerHTML;
            y.innerHTML = val;
        }
    }
}

function reset_X(){
    // Reset the params shown in the X_n box.
    document.getElementById('X_el').innerHTML = 'Mean'
    document.getElementById('N').innerHTML = '';
    document.getElementById('X_n').style.backgroundColor = "#ffffff";
    for (key in keys){
        document.getElementById("X_" + keys[key]).innerHTML = document.getElementById("av_" + keys[key]).innerHTML
    }
}

function set_av(val){
    document.getElementById('av_el').innerHTML = 'Mean'
    for (key in keys){
        document.getElementById("av_"+keys[key]).innerHTML = val[key]
    }
}

function show_X(){
    // Show the X_n box.
    document.getElementById('X_n').style.visibility = "visible";
}

function hide_X(){
    // Hide the X_n box.
    document.getElementById('X_n').style.visibility = "hidden";
}

function humanize(size) {
	var units = ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB'];
	var ord = Math.floor(Math.log(size) / Math.log(1024));
	ord = Math.min(Math.max(0, ord), units.length - 1);
	var s = Math.round((size / Math.pow(1024, ord)) * 100) / 100;
	return s + ' ' + units[ord];
}

function dojoTour_guidedtour() {
    var intro = introJs();
    intro.setOptions({
      steps: [
        {
          intro: "Welcome to the PseudoDojo! Let me explain how to use the website."
        },
        {
          element: '#TYP',
          intro: 'Here you select the type of pseudopotential. SR stands for scalar relativistic, FR for fully relativistic (including SOC). '+
                 'The options for xc, accuracy and format are adjusted based on your choice.'
        },
        {
          element: '#XCF',
          intro:  "In this selector you can pick one of the available exchange-correlation functionals. " +
                  "Have a look at the F.A.Q. if your fuctional of choice is not there."
        },
        {
          element: '#ACC',
          intro:  "In this selector you can select one of the available accuracy levels. " +
                  "Have a look at the F.A.Q. for a detailed description."
        },
        {
          element: '#FMT',
          intro:  "In this selector you can pick the format of the pseudopotential file. " +
                  "PSP8 for ABINIT, UPF2 for Quantum Espresso, PSML1.1 is supported by Siesta. " +
                  "When you select HTML, clicking the elements will display a full report of all the tests we performed. " +
                  "Finally djrepo will give you all the numerical results of the tests in json format."
        },
        {
          element: '#X_n',
          intro:  "As long as you don't hover one of the elements, this box shows the average values for the table you selected. " +
                  "Once you hover the elements, it shows the values for that element. "
        },
        {
          element: "#X_hl",
          intro:  "This is the low cutoff energy hint (Ha). Good for a quick calculation or as a starting point for your convergence studies."
        },
        {
          element: "#X_hn",
          intro:  "This is the normal cutoff energy hint (Ha). A good guess for high-throughput calculations."
        },
        {
          element: "#X_hh",
          intro:  "This is the high cutoff energy hint (Ha), beyond this value you should not observe significant changes in the results anymore."
        },
        {
          element: "#X_nv",
          intro:  "The number of valence shells."
        },
        {
          element: "#X_d",
          intro:  "The results of the delta gauge test (meV). This is the integral between the equation of state calculated using the pseudo potential and a reference all electron equation of state."
        },
        {
          element: "#X_dp",
          intro:  "The normalized delta gauge."
        },
        {
          element: "#X_gb",
          intro:  "The gbrv fcc and bcc average (%). This is the relative error in the lattice parameter with respect to reference all electrons results."
        },
        {
          element: "#silicon",
          intro:  "You can now click all the elements in the table to download or view the selected file for a single element. If the box turns green the file is available, if it turns red we are still working on it."
        },
        {
          element: ".download_button",
          intro:  "Alternatively, with the download button you can get a tar of the full table, always one pseudopotential per element."
        },
        {
          element: "#papers",
          intro:  "A list of papers using PseudoDojo pseudopotentials. Did you use them? Send us the DOI and we'll add yours as well."
        },
        {
          element: ".logo",
          intro:  "Finally, if you want to learn the periodic table by heart try clicking here. (p.s. Don't try to download Oganesson, bad things may happen.)"
        }
      ],
      showProgress: true,
      overlayOpacity: 0.3
    });

    //var remove_glow = function() {
    //    console.log('exiting');
    //    $("#guided-tour-button").removeClass("glow");
    //};
    //intro.onexit(remove_glow);
    //intro.oncomplete(remove_glow);

    intro.start();
}

function dynamicdropdown(listindex){
  // Set the values of the XC/Accuracy/Format widgets given the pseudo type.
  console.log('dynamic dropdown: setting', listindex)
  document.getElementById("ACC").length = 0;
  document.getElementById("XCF").length = 0;
  document.getElementById("FMT").length = 0;

  switch (listindex)
  {
    case "paw" :
      document.getElementById('warning_box').innerHTML = "";
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      document.getElementById("ACC").options[1]=new Option("stringent","stringent");
      document.getElementById("XCF").options[1]=new Option("LDA","pw");
      document.getElementById("XCF").options[0]=new Option("PBE","pbe");
      document.getElementById("FMT").options[0]=new Option("xml","xml");
      break;

    case "nc-sr" :
      document.getElementById('warning_box').innerHTML = "";
      set_warning(' this version is outdated')
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      document.getElementById("ACC").options[1]=new Option("stringent","stringent");
      document.getElementById("XCF").options[2]=new Option("LDA","pw");
      document.getElementById("XCF").options[0]=new Option("PBE","pbe");
      document.getElementById("XCF").options[1]=new Option("PBEsol","pbesol");
      document.getElementById("FMT").options[0]=new Option("psp8","psp8");
      document.getElementById("FMT").options[1]=new Option("upf","upf");
      document.getElementById("FMT").options[2]=new Option("html","html");
      document.getElementById("FMT").options[3]=new Option("djrepo","djrepo");
      break;

    case "nc-sr-04" :
      document.getElementById('warning_box').innerHTML = "";
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      document.getElementById("ACC").options[1]=new Option("stringent","stringent");
      document.getElementById("XCF").options[2]=new Option("LDA","pw");
      document.getElementById("XCF").options[0]=new Option("PBE","pbe");
      document.getElementById("XCF").options[1]=new Option("PBEsol","pbesol");
      document.getElementById("FMT").options[0]=new Option("psp8","psp8");
      document.getElementById("FMT").options[1]=new Option("upf","upf");
      document.getElementById("FMT").options[2]=new Option("psml","psml");
      document.getElementById("FMT").options[3]=new Option("html","html");
      document.getElementById("FMT").options[4]=new Option("djrepo","djrepo");
      break;

    case "nc-fr-04" :
      document.getElementById('warning_box').innerHTML = "";
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      document.getElementById("ACC").options[1]=new Option("stringent","stringent");
      //document.getElementById("XCF").options[2]=new Option("LDA","pw");
      document.getElementById("XCF").options[0]=new Option("PBE","pbe");
      document.getElementById("XCF").options[1]=new Option("PBEsol","pbesol");
      document.getElementById("FMT").options[0]=new Option("psp8","psp8");
      document.getElementById("FMT").options[1]=new Option("upf","upf");
      document.getElementById("FMT").options[2]=new Option("psml","psml");
      document.getElementById("FMT").options[3]=new Option("html","html");
      document.getElementById("FMT").options[4]=new Option("djrepo","djrepo");
      break;

    case "nc-sr-04-3plus" :
      set_warning(" this table contains Lanthanide potentials for use in the 3+ configuration only. <b>They all have the f-electrons frozen in the core.</b> The hints are based on the convergence of the nitride lattice parameter, see the report under format:html for details.");
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      // document.getElementById("ACC").options[1]=new Option("stringent","stringent");
      // document.getElementById("XCF").options[]=new Option("LDA","pw");
      document.getElementById("XCF").options[0]=new Option("PBE","pbe");
      // document.getElementById("XCF").options[]=new Option("PBEsol","pbesol");
      document.getElementById("FMT").options[0]=new Option("psp8","psp8");
      document.getElementById("FMT").options[1]=new Option("upf","upf");
      document.getElementById("FMT").options[2]=new Option("psml","psml");
      document.getElementById("FMT").options[3]=new Option("html","html");
      document.getElementById("FMT").options[4]=new Option("djrepo","djrepo");
      break;

    case "core" :
      // TODO or perhaps add new format and handle file download.
      document.getElementById('warning_box').innerHTML = "";
      document.getElementById("ACC").options[0]=new Option("","standard");
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      document.getElementById("ACC").options[0]=new Option("standard","standard");
      document.getElementById("XCF").options[0]=new Option("PBE","pbe");
      document.getElementById("FMT").options[2]=new Option("FC","fc");
  }
  return true;
}

function chaos() {
    localStorage.setItem('animate', 1)
    $('.plugin').removeClass('anim');
    $('.plugin').removeClass('chaos');
    setTimeout("$('.plugin').addClass('chaos')",10)
    var plugins = document.querySelectorAll(".plugin");
    for (var i = 0; i < 118; i++) {
      var plugin = plugins[i];
      animatePlugin(plugin);
    }
    function animatePlugin(plugin) {
      var xMax = 500;
      var yMax = 500;
      var x1 = Math.random() - 0.5;
      x1 = x1 * xMax;
      var x2 = Math.random() - 0.5;
      x2 = x2 * xMax;
      var y1 = Math.random() - 0.5;
      y1 = y1 * yMax;
      var y2 = Math.random() - 0.5;
      y2 = y2 * yMax;

      plugin.keyframes = [{
        opacity: 1,
        transform: "translate3d(" + x1 + "px, " + y1 + "px, 0px)"
      }, {
        opacity: 0.2,
        transform: "translate3d(" + x2 + "px, " + y2 + "px, 0px)"
      }, {
        opacity: 0.2,
        transform: "translate3d(" + -x1 + "px, " + -y1 + "px, 0px)"
      }, {
        opacity: 1,
        transform: "translate3d(" + -x2 + "px, " + -y2 + "px, 0px)"
      }];

      plugin.animProps = {
        duration: 2000 + Math.random() * 4000,
        easing: "ease-out",
        iterations: 1
      }
    var animationPlayer = plugin.animate(plugin.keyframes, plugin.animProps);
    }
}
