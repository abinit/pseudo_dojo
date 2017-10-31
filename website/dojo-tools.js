/*
function.prototype.method = function (name func) {
    this.prototype[name] = func
    return this
}
*/

var keys = ['hh', 'hl', 'hn', 'nv', 'd', 'dp', 'gb'];
var els = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne','Na', 'Mg', "Al", "Si", 'P', 'S', 'Cl', 'Ar','K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'];

function set_info(info) {
    var averages = [0,0,0,0,0,0,0];
    var sums = [0,0,0,0,0,0,0];
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
            else{
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
    if (sums[0] > 0){
            //alert(summary);
            set_av(averages);
            reset_X()
            }
}

function loadJSON(file, callback) {
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

function load_set_info() {
    acc = document.getElementById('ACC').value;
    xcf = document.getElementById('XCF').value;
    type = document.getElementById('TYP').value;
    set_info({});
    var file = type + '_' + xcf + '_' + acc + '.json';
    loadJSON(file, function(response) {
    var info = JSON.parse(response);
    set_info(info);
    });
}

function set_X(elm){
    if (els.indexOf(elm)>=0){
        var x = document.getElementById('X_el');
        x.innerHTML = elm;
        for (key in keys){
            var id_key = 'X_' + keys[key];
            var id_key_in = elm + '_' + keys[key];
            var x = document.getElementById(id_key_in);
            var y = document.getElementById(id_key);
            var val = x.innerHTML
            y.innerHTML = val;
        }
    }
}

function reset_X(){
    document.getElementById('X_el').innerHTML = 'Mean'
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
    document.getElementById('X_n').style.visibility="visible";
}

function hide_X(){
    document.getElementById('X_n').style.visibility="hidden";
}


var x = 5;

var p_text = '';
p_text += 'hello world first';

function to_red() {
    var x = document.getElementById("par-1");
    x.style.color = "red";
    x.innerHTML = "test changed to red"
}

function to_black() {
    var x = document.getElementById("par-1");
    x.style.color = "black";
    x.innerHTML = "test changed to black"
}

function to_color() {
    var color = document.getElementById("color_in").value
    //alert('setting paragraph to ' + color);
    var x = document.getElementById("par-1");
    x.style.color = color;
    x.innerHTML = "test changed to custom color"
}

function to_color_arg(cc) {
    //alert('setting paragraph to ' + cc);
    var x = document.getElementById("par-2");
    x.style.color = cc;
    x.innerHTML = "test changed to custom color via argument"
}

function product(x,y) {
    var p = x * y;
    var r = document.getElementById('result-1');
    r.innerHTML = 'the product of ' + x + ' and ' + y + ' is ' + p
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
          intro: 'Here you select the type of pseudopotential. '+
                 'The options for xc, accuracy and format are adjusted based on your choice here.'
        },
        {
          element: '#XCF',
          intro:  "In this selector you can pick one of the available exchange correlation functionals. " +
                  "Have a look at the F.A.Q. if your fuctional of choice is not there."
        },
        {
          element: '#ACC',
          intro:  "In this selector you can pick one of the accuracies. " +
                  "Have a look at the F.A.Q. for a detailed description."
        },
        {
          element: '#FMT',
          intro:  "In this selector you can pick the format of the pseudopotential file. " +
                  "PSP8 for ABINIT, UPF (UPF2) for quantum espresso. " +
                  "HTML will give you a full report on the tests. " +
                  "djrepo will give you the full numerical results of the tests."
        },
        {
          element: '#X_n',
          intro:  "As long as you don't hover one of the elements this box shows the average values for the table you selected. " +
                  "Once you hover the elements it shows the values for that element. "
        },
        {
          element: "#X_hl",
          intro:  "low cutoff energy hint (Ha)"
        },
        {
          element: "#X_hn",
          intro:  "normal cutoff energy hint (Ha)"
        },
        {
          element: "#X_hh",
          intro:  "high cutoff energy hint (Ha)"
        },
        {
          element: "#X_nv",
          intro:  "number of valence shells"
        },
        {
          element: "#X_d",
          intro:  "delta gauge (meV):&#013;Integral between the equation of state calculated using the pseudo potential &#013;and a reference all electron equation of state."
        },
        {
          element: "#X_dp",
          intro:  "normalized delta gauge:&#013;Normalized version of the delta gauge."
        },
        {
          element: "#X_gb",
          intro:  "gbrv fcc bcc average (%):&#013;Relative error in the lattice parameter with respect to reference."
        },
        {
          element: "#silicon",
          intro:  "You can now click all the elements in the table to download single elements. If the box turns green the file is available, if it turns red.... "
        },
        {
          element: "#download_button",
          intro:  "Alternatively, with the download button you can now get a tar of the full table, always one pseudopotential per element."
        },
        {
          element: ".logo",
          intro:  "Finally, if you want to learn the periodic table by hard, click here."
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