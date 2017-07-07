/*
function.prototype.method = function (name func) {
    this.prototype[name] = func
    return this
}
*/

var keys = ['hh', 'hl', 'hn', 'nv', 'd', 'dp', 'gb'];
var els = ['H','He', 'Li', 'Be', 'B'];

function set_info(info) {
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
            x.innerHTML = val;
        }
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
    type = document.getElementById('REL').value;
    set_info({});
    var file = type + '_' + xcf + '_' + acc + '.json';
    alert(file)
    loadJSON(file, function(response) {
    var info = JSON.parse(response);
    set_info(info);
    });
}

function set_X(elm){
    if (elm != 'X'){
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
    document.getElementById('X_el').innerHTML = 'X'
    document.getElementById("X_hl").innerHTML = 'low'
    document.getElementById("X_hn").innerHTML = 'normal'
    document.getElementById("X_hh").innerHTML = 'high'
    document.getElementById("X_nv").innerHTML = 'nv'
    document.getElementById("X_d").innerHTML = 'delta'
    document.getElementById("X_dp").innerHTML = 'delta1'
    document.getElementById("X_gb").innerHTML = 'gbrv'
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
