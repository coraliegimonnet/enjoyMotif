/* Absolute pathway of position file */
var path = window.location.pathname.split('/'),
parent = path[path.length-2];

path.splice(0,1);

var pathway = "";

for (var i = 0 ; i < path.length - 1 ; i ++)
{
    pathway += path[i]+"/";
}


/* Transcript's name which be drawing */

var gene = prompt("Gene Symbol : \n(ex : RpL36)");
if (/\S+\(\S+\)(\S+)/.exec(gene))
{
    //console.log(RegExp.$1);
    var regexend = RegExp.$1;
    var regex_parenthesis = new RegExp (regexend, 'i');
}
else
{
    var regex = new RegExp ("^"+gene, 'i');
}


/* All functions used for drawing and parsing file input */

function getXMLHttpRequest()
{
    var xmlhttp = null;
    if (window.XMLHttpRequest || window.ActiveXObject) {
        if (window.ActiveXObject) {
            try {
                xmlhttp = new ActiveXObject("Msxml2.XMLHTTP");
            }
            catch(e) {
                xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
            }
        }
        else {
            xmlhttp = new XMLHttpRequest();
        }
    }
    else {
        alert("Your browser does not support the XMLHTTPRequest object.");
        return null;
    }
    return xmlhttp;
}

function getRSS(file)
{
    var rss = getXMLHttpRequest();
    rss.open("GET", file , false);
    rss.send(null);
    var line = rss.responseText.split(/\n/g);
    return(line);
}

/* Object creation */

function myRBP(table)
{
    var RBP = {
        utr : table[0],
        long : table[3],
        name : table[1],
        motif : table[2],
        start : table[4],
        end : table[5]
    };
    return RBP;
}

function MyPolyA(table)
{
    var polyA = {
        utr : table[0],
        long : table[3],
        name : table[1],
        motif : table[2],
        start : table[4],
        end : table[5]
    };
    return polyA;
}

function cleanArray(array)
{
    var i, j, len = array.length, out = [], obj = {};
    for (i = 0; i < len; i++) {
        obj[array[i]] = 0;
    }
    for (j in obj) {
        out.push(j);
    }
    return out;
}

function ProtBinding (prot,col)
{
    var RBP = {
        name : prot,
        color : col
    };
    return RBP;
}

function MyMIR (table)
{
    var mir = {
        utr : table[0],
        name : table[1],
        start : table[2],
        end : table[3],
        sequence : table[4]
    };
    return mir;
}

function MyTranscript (table)
{
    var transcript = {
        utr : table[0],
        long : table[3],
        onebp : 1000 / table[3]
    };
    return transcript;
}

/* Remove duplicate in Object Array */
function removeDuplicates(originalArray, objKey) 
{
    var trimmedArray = [];
    var values = [];
    var value;
    for(var i = 0; i < originalArray.length; i++) 
    {
        value = originalArray[i][objKey];

        if(values.indexOf(value) === -1) 
        {
            trimmedArray.push(originalArray[i]);
            values.push(value);
        }
    }
    return trimmedArray;
}

/* miRNAs */

var tableMIR = getRSS("file:///"+pathway+"mirna_alignment.txt");
var outputMIR = [];

for (var i = 0 ; i < tableMIR.length ; i++)
{
    if (typeof regex_parenthesis !== "undefined")
    {
        
        if (regex_parenthesis.exec(tableMIR[i]))
        {
            outputMIR.push(tableMIR[i]);
            
        }
    }
    else if (tableMIR[i].match(regex))
    {
        outputMIR.push(tableMIR[i]);
        //console.log("blabla");
    }
    
}

var miRNAs = new Array;

for (var i = 0 ; i < outputMIR.length ; i++)
{
    var drawMIR = outputMIR[i].split(",");
    var mir = new MyMIR(drawMIR);
    miRNAs.push(mir);
    //console.log(mir);
}



/* 3'UTR and site */

var table = getRSS("file:///"+pathway+"position.csv");
var output=[];

for (var i = 0 ; i < table.length ; i++)
{
    if (typeof regex_parenthesis !== "undefined")
    {
        if (regex_parenthesis.exec(table[i]))
        {
            output.push(table[i]);
        }
    }
    else if (table[i].match(regex))
    {
        output.push(table[i]);
    }
}

var TableTranscript = new Array;
var PolyA = new Array;
var TableRBP = new Array;
var protein = new Array;
var length;


for (var i = 0 ; i < output.length ; i++)
{
    var draw = output[i].split(",");
    var transcript = new MyTranscript(draw);
    TableTranscript.push(transcript);
    if (/poly\(A\)/i.test(draw[1]))
    {
        var polyA = new MyPolyA(draw);
        PolyA.push(polyA);
        //console.log(polyA);
    }
    else
    {
        protein.push(draw[1]);
        var RBP = new myRBP(draw);
        TableRBP.push(RBP);
        //console.log(RBP);
    }
}

var TrimmedTranscript = removeDuplicates (TableTranscript,'utr');
// for (var i = 0 ; i < TableTranscript.length ; i++ )
// {
//     console.log(TableTranscript[i]);
// }

var sortProtein = cleanArray(protein);
var ProtColor = new Array;

for (var i = 0 ; i < sortProtein.length ; i++)
{
    var RBPcolor = prompt("Which color for this binding protein : " + sortProtein[i] +"? \n(Warning : yellow and orange are used for PAS)");
    var protColor = new ProtBinding(sortProtein[i],RBPcolor);
    ProtColor.push(protColor);
}


for (var j = 0 ; j < TrimmedTranscript.length ; j++)
{
    for (var i = 0 ; i < PolyA.length ; i++)
    {
        if (TrimmedTranscript[j].utr === PolyA[i].utr)
        {
            var x = PolyA[i].start * TrimmedTranscript[j].onebp;
            PolyA[i].startpx = Math.round(x);
            var y = PolyA[i].end * TrimmedTranscript[j].onebp;
            PolyA[i].endpx = Math.round(y);
            //console.log(PolyA[i]);
        }
        
    }
    for (var i = 0 ; i < TableRBP.length ; i++)
    {
        for (var l = 0 ; l < ProtColor.length ; l++)
        {
            if (ProtColor[l].name === TableRBP[i].name)
            {
                var c = ProtColor[l].color;
                TableRBP[i].color = c;
            }
        }
        if (TrimmedTranscript[j].utr === TableRBP[i].utr)
        {
            var x = TableRBP[i].start * TrimmedTranscript[j].onebp;
            TableRBP[i].startpx = Math.round(x);
            var y = TableRBP[i].end * TrimmedTranscript[j].onebp;
            TableRBP[i].endpx = Math.round(y);
            //console.log(TableRBP[i]);
        }
    }

    var StartY = 390;
    for (var i = 0 ; i < miRNAs.length ; i++)
    {
        if (TrimmedTranscript[j].utr === miRNAs[i].utr)
        {
            var x = miRNAs[i].start * TrimmedTranscript[j].onebp;
            miRNAs[i].startpx = Math.round(x);
            var y = miRNAs[i].end * TrimmedTranscript[j].onebp;
            miRNAs[i].endpx = Math.round(y);
            miRNAs[i].starty = StartY;
            StartY+=20;
            //console.log(miRNAs[i]);
        }
    }
}


/* Canvas design */

window.onload = function()
{
    var canvas = document.getElementById('mon_canvas');
    if(!canvas)
    {
        alert("Unable to retrieve the canvas");
        return;
    }

    var context = canvas.getContext('2d');
    if(!context)
    {
        alert("Unable to retrieve the context of canvas");
        return;
    }
    context.font = "30px Helvetica";
    context.fillText(gene,500,40);
    var A = 0;
    for (var x = 0 ; x < TrimmedTranscript.length ; x++)
    {
        context.globalAlpha = 1;
        //Title
        context.fillStyle = "black";
        context.font="15px Helvetica";
        context.fillText(TrimmedTranscript[x].utr+" (" + TrimmedTranscript[x].long + " bp)" , 5 , 300 + A);
        //UTR
        context.fillStyle = "rgba(117,177,255,1)";
        context.fillRect(10,320+A,1000,30);
        context.fillStyle = "black";
        context.font = "15px Helvetica";
        context.fillText("3'UTR",2,365+A);
        //Caption RBP
        var y = 110;
        var text = 122;
        for (var i = 0 ; i < ProtColor.length ; i++)
        {
            context.fillStyle = ProtColor[i].color;
            context.fillRect(10,y,15,15);
            context.strokeRect(10,y,15,15);
            context.fillStyle = "black";
            context.font = "15px Helvetica";
            context.fillText(ProtColor[i].name,40,text);
            y+=30;
            text+=30;
        }
        for (var i = 0 ; i < PolyA.length ; i++)
        {
            if (/non-canonical/.test(PolyA[i].name))
            {
                    //Caption
                context.fillStyle = "rgba(233,174,0,1)";
                context.fillRect(10,80,15,15);
                context.strokeRect(10,80,15,15);
                context.fillStyle = "black";
                context.font = "15px Helvetica";
                context.fillText("Polyadenylation Signal (non-canonical)",40,92);

                context.fillStyle = "rgba(233,174,0,0.7)";
            }
            else
            {
                //Caption
                context.fillStyle = "rgba(255,255,0,1)";
                context.fillRect(10,50,15,15);
                context.strokeRect(10,50,15,15);
                context.fillStyle = "black";
                context.font = "15px Helvetica";
                context.fillText("Polyadenylation Signal (canonical)",40,62);

                context.fillStyle = "rgba(255,255,0,0.7)";
            }
            if (PolyA[i].utr === TrimmedTranscript[x].utr)
            {
                canvas.addEventListener("click", click_manage_polyA, false);
                var startX = 10 + PolyA[i].startpx;
                var endX = PolyA[i].endpx - PolyA[i].startpx;
                context.globalCompositeOperation = "source-over";
                context.fillRect(startX,320+A,endX,30);
            }

        }
        //RBP
        for (var i = 0 ; i < TableRBP.length ; i++)
        {
            if (TableRBP[i].utr === TrimmedTranscript[x].utr)
            {
                canvas.addEventListener("click", click_manage_RBP, false);
                var startX = 10 + TableRBP[i].startpx;
                var endX = TableRBP[i].endpx - TableRBP[i].startpx;
                //context.fillStyle = "rgba(215,44,44,0.7)";
                context.globalAlpha = 0.7;
                context.fillStyle = TableRBP[i].color;
                context.globalCompositeOperation = "source-over";
                context.fillRect(startX,320+A,endX,30);
            }
            
        }
        //miRNA
        context.globalAlpha = 1;
        context.fillStyle = "black";
        context.font = "15px Helvetica";
        context.fillText("miRNAs",2,390+A);
        var startRect = 35;
        for (var i = 0 ; i < miRNAs.length ; i++)
        {
            if (TrimmedTranscript[x].utr === miRNAs[i].utr)
            {
                var startY = miRNAs[i].starty;
                canvas.addEventListener("click", click_manage_mir, false);
                var startX = 10 + miRNAs[i].startpx;
                var endX = miRNAs[i].endpx - miRNAs[i].startpx;
                context.fillStyle = "rgba(0,0,0,1)";
                context.globalCompositeOperation = "source-over";
                context.fillRect(startX,startY+A,endX,10);
                startRect += 20;
            }
            
        }
        context.fillStyle = "black";
        context.strokeRect(0,375+A,1000,startRect);
        A+=375;
    }
    
    //Save canvas as image
    function dlCanvas()
    {
        var dt = canvas.toDataURL('image/png');
        dt = dt.replace(/^data:image\/[^;]*/, 'data:application/octet-stream');
        dt = dt.replace(/^data:application\/octet-stream/, 'data:application/octet-stream;headers=Content-Disposition%3A%20attachment%3B%20filename=Canvas.png');

        this.href = dt;
    };
    document.getElementById("dl").addEventListener('click', dlCanvas, false);
}

function click_manage_polyA (e)
{
    var i, pol;
    var box = this.getBoundingClientRect(); // "this" est le canvas
    var x_clic = e.clientX - box.left;
    var y_clic = e.clientY - box.top;
    var A = 0;
    for (var j = 0 ; j < TrimmedTranscript.length ; j++)
    {
        for (i = 0 ; i < PolyA.length ; i++)
        {
            if (PolyA[i].utr === TrimmedTranscript[j].utr)
            {
                pol = PolyA[ i ];
                if ((x_clic >= 10 + pol.startpx ) && (x_clic <= 10 + pol.endpx) && (y_clic >= 320 + A) && (y_clic <= 350 + A))
                {
                    alert(pol.name+" : " + pol.motif+"\nPosition : "+pol.start+"-"+pol.end+"\n");
                }
            }
        }
        A+=375;
    }
}

function click_manage_RBP (e)
{
    var i, rbp;
    var box = this.getBoundingClientRect(); // "this" est le canvas
    var x_clic = e.clientX - box.left;
    var y_clic = e.clientY - box.top;
    var A = 0;
    for (var j = 0 ; j < TrimmedTranscript.length ; j++)
    {
        for (i = 0 ; i < TableRBP.length ; i++)
        {
            if (TableRBP[i].utr === TrimmedTranscript[j].utr)
            {
                rbp = TableRBP[ i ];
                if ((x_clic >= 10 + rbp.startpx ) && (x_clic <= 10 + rbp.endpx) && (y_clic >= 320+A) && (y_clic <= 350+A))
                {
                    alert(rbp.name+" : " + rbp.motif+"\nPosition : "+rbp.start+"-"+rbp.end+"\n");
                }
            }
        }
        A += 375;
    }
}

function click_manage_mir (e)
{
    var i, miR;
    var box = this.getBoundingClientRect(); // "this" est le canvas
    var x_clic = e.clientX - box.left;
    var y_clic = e.clientY - box.top;
    var A = 0;
    for (var j = 0 ; j < TrimmedTranscript.length ; j++)
    {
        for (i = 0 ; i < miRNAs.length ; i++)
        {
            if (miRNAs[i].utr === TrimmedTranscript[j].utr)
            {
                miR = miRNAs[ i ];
                if ((x_clic >= 10 + miR.startpx ) && (x_clic <= 10 + miR.endpx) && (y_clic >= miR.starty + A) && (y_clic <= miR.starty+10+A))
                {
                    alert(miR.name+" : " + miR.sequence+"\nPosition : "+miR.start+"-"+miR.end+"\n");
                }
            }
        }
        A += 375;
    }
}



