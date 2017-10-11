/*
function.prototype.method = function (name func) {
    this.prototype[name] = func
    return this
}
*/

var keys = ['hh', 'hl', 'hn', 'nv', 'd', 'dp', 'gb'];
var els = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne','Na', 'Mg', "Al", "Si", 'P', 'S', 'Cl', 'Ar','K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'];

function set_info(info, animate) {
    var averages = [0,0,0,0,0,0,0];
    var sums = [0,0,0,0,0,0,0];
    if (animate === 1){
        console.log('added animating');
        $('.plugin').removeClass('anim');
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

function store_available_files() {
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
    if (els.indexOf(elm)>=0){
        document.getElementById('N').innerHTML = n;
        var x = document.getElementById('X_n');
        x.style.backgroundColor = color;
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
                  "When you select HTML clicking the elements will display a full report on all the tests we performed. " +
                  "Finally djrepo will give you the full numerical results, in json format, of the tests."
        },
        {
          element: '#X_n',
          intro:  "As long as you don't hover one of the elements this box shows the average values for the table you selected. " +
                  "Once you hover the elements it shows the values for that element. "
        },
        {
          element: "#X_hl",
          intro:  "This is the low cutoff energy hint (Ha) good for a quick first calculation or as a starting point for your convergence studies."
        },
        {
          element: "#X_hn",
          intro:  "This is the normal cutoff energy hint (Ha) a good guess for high through put calculations."
        },
        {
          element: "#X_hh",
          intro:  "This is the high cutoff energy hint (Ha), beyond this value you should not observe any real changes to the results anymore."
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
          intro:  "Alternatively, with the download button you can now get a tar of the full table, always one pseudopotential per element."
        },
        {
          element: ".logo",
          intro:  "Finally, if you want to learn the periodic table by hard try clicking here."
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