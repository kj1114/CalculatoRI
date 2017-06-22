/* Created by Krstina Janjic on 20.5.17 */

//Main JavaScript File

//function triggered on the click of the form submit

function submitFunction(){

    //REMOVES CHILD NODE OF THE CHART ELEMENT WHEN FORM IS RE-SUBMITTED
    //THIS ENABLES THE USER TO INPUT NEW DATA WITHOUT HAVING TO REFRESH THE PAGE

	var parent = document.getElementById('insertChart');
	var child = document.getElementById('morris-line-chart');
	parent.removeChild(child);

	//FETCHES THE CONTENT OF INPUT FIELDS
    //INSERTS IT INTO -formData- OBJECT

    var conc, ref, rad, lmd, f;
    var formData = new FormData();

    //in case an input file is not uploaded, oxygenated and deoxygenated options needed
    if (document.getElementById('radio1').checked){   //oxygenated button is selected
        rad = true;
    }
    else{                                                   //deoxygenated button is selected
        rad = false;
    }
    conc = document.getElementById('concentration').value;
    ref = document.getElementById('Nreference').value;
    lmd = document.getElementById('lambdaReference').value;
    f = document.getElementById('file1').files[0];

    //insert the input fields values into formData object
    formData.append('concentration', conc);
    formData.append('Nreference', ref);
    formData.append('radio1', rad);
    formData.append('lambdaReference', lmd);
    formData.append('file1', f);

    //Test if 'Oxygenated' button is checked

    console.log (document.getElementById('radio1').checked);

    //SEND THE INPUT DATA OBJECT TO THE SERVER
    //makes use of XMLHttpRequest built-in object to interact with the server

    //tell the server that you are sending data, hence POST
    //specify URL of the server
    var req = new XMLHttpRequest();

    req.open('POST', 'http://localhost:8080/callCPPFunc');
    req.send(formData);

    req.onreadystatechange = function() {

        if (this.readyState == 4 && this.status == 200) {
            var recText = JSON.parse(this.responseText);
            var chart = document.createElement('div');
            chart.setAttribute('id', 'morris-line-chart');
            parent.appendChild(chart);
            Morris.Line({
		        element: 'morris-line-chart',
		        data: recText,
		        xkey: 'wavelength',
		        ykeys: ['Re_n'],
		        ymax: 'auto',
		        ymin: 'auto',
                verticalGrid: true,
		        labels: ['Re_n'],
                lineWidth: 1.5,
		        pointSize: 0.05,
                yLabelFormat: function (y) {
		            return y.toPrecision(4).toString();
		            },
                xLabelFormat: function (x){
		            return x.src.wavelength.toPrecision(4).toString();
		            },
                hoverCallback: function (index, options, content, recText) {
		            var verticalAxis = 'Real Part RI: ' + recText.Re_n.toPrecision(4);
		            var horizontalAxis = 'Wavelength: ' + recText.wavelength.toPrecision(4) + ' nm';
                    return [verticalAxis, horizontalAxis].join('<br/>');
                },
                parseTime: false,
		        resize: true
		    });
            var getRe_n = [];
            var getWavelength = [];

            for (i = 0; i < recText.length; i++) {
                getRe_n.push(recText[i].Re_n.toPrecision(8));
                getWavelength.push(recText[i].wavelength.toPrecision(4))
            }

            var contentTxt = "";
            for (i = 0; i < recText.length; i++) {
                contentTxt+= getWavelength[i] +  "   " + getRe_n[i] + "\n";
            }
            var maxes = [];
            var condition = [];
            var index_find = [];

            function checkTrue(age) {
                return age = true;
            }

            for (var i = 1; i < getRe_n.length - 1; ++i) {
                for (var j = 1; j < 2; ++j) {
                    condition.push(getRe_n[i - j] < getRe_n[i] && getRe_n[i] > getRe_n[i + j]);
                }
                if (condition.every(checkTrue)){
                    index_find.push(i);
                }
                condition = [];
            }
            for (var m =0; m < index_find.length; ++m){
                maxes.push(getRe_n[index_find[m]]);
            }

            var dataStrTxt = "data:text/txt;charset=utf-8," + encodeURIComponent(contentTxt);
            var dlAnchorElemTxt = document.getElementById('export-text-file');
            dlAnchorElemTxt.setAttribute("href",     dataStrTxt     );
            dlAnchorElemTxt.setAttribute("download", "RealPartRefractiveIndexVsLambda.txt");
        }

        if (this.readyState == 4 && this.status == 500) {
            var chart = document.createElement('div');
            chart.setAttribute('id', 'morris-line-chart');
            parent.appendChild(chart);
        }
    }
}

$('input[type=file]').change(function(){
    if(document.getElementById("file1").files.length === 0){
        $("#buttons-section").show();
        $("#concentration").show();
    }
    else {
        $("#buttons-section").hide();
        $("#concentration").hide();
    }
});
