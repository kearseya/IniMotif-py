from KmerKounter import identifier
from KmerKounter import numofruns
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import filenames1

from PositionBias import TSeqNums
from PositionBias import LSeqNums
from PositionBias import Barcodevalues
from PositionBias import numofkmers
from PositionBias import numofuniquekmers
from PositionBias import numoftfbs
from PositionBias import numoftfbsseq
from PositionBias import seqbias


def setidentity():
    global htmlname
    htmlname = identifier+'.html'
setidentity()
#print(identifier, htmlname)



Html_file=open(str(htmlname), "w")

title = """
<!DOCTYPE html>
<html lang="en">
<head>
<title>CSS Template</title>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
* {
  box-sizing: border-box;
}

body {
  font-family: Arial, Helvetica, sans-serif;
}

p {
  font-size: 100px
}

h1 {
  font-size: 200px
}

h2 {
  font-size: 35px
}


header {
  background-color: #666;
  padding: 10px;
  text-align: center;
  font-size: 35px;
  color: white;
}
</style>
</head>
<header>
<p> IniMotif </p>
<h2>A pipeline for motif discovery!</h2>
</header>
<br>
"""

Html_file.write(title)

def runheader():
    runheaders = []
    for x in range(0, numofruns+1):
        string="<h1>Run "+str(x)+"</h1>"
        runheaders.append(string)
    return runheaders

runheaders = runheader()


def results():
    html_strs = []
    for z in range(0, numofruns+1):
        html_strs.append({})
        for i in range(mink, maxk+1):
            html_strs[z].update({i:[]})

    for x in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            string = """<tr><td><p style="background-color: lightgrey; color: #333; padding: 30px;">K = """+str(k)+"""</p><img src="figures/logo_"""+str(identifier)+"_"+str(x)+"_"+str(k)+""".png"><img src="figures/hamdist_"""+str(identifier)+"_"+str(x)+"_"+str(k)+""".png" width=2000px hight=1200px><img src="figures/pos_"""+str(identifier)+"_"+str(x)+"_"+str(k)+""".png" width=2000px hight=1200px></td></tr>"""
            html_strs[x][k].append(string)

    return html_strs

html_strs = results()


def kmerfrequency():
    html_kmerfreq = {}
    for k in range(mink, maxk+1):
        html_kmerfreq.update({k:[]})
        #print(html_kmerfreq)
        #html_kmerfreq.update({k:[]})
        string = """<p style="background-color: lightgrey; color: #333; padding: 30px;">K = """+str(k)+"""</p>"""+"""<img src="figures/kmerfreq_"""+str(identifier)+"_"+str(k)+""".png" width=3000px height=3000px>"""
        html_kmerfreq[k].append(string)

    return html_kmerfreq

html_kmerfreq = kmerfrequency()
#print(html_kmerfreq)


def formatter():
    for r in range(1,numofruns+1):
        Html_file.write(runheaders[r])
        for k in range(mink, maxk+1):
            Html_file.write(html_strs[r][k][0])

def kmerfreqformatter():
    Html_file.write("""<h1 style="background-color: #3c3c3c; padding: 20px; color: white;"> Kmer frequency </h1>""")
    for i in range(mink, maxk+1):
        Html_file.write(html_kmerfreq[i][0])

formatter()
kmerfreqformatter()

def tablemaker():
    Html_file.write("""<h1 style="background-color: #3c3c3c; padding: 20px; color: white;"> Numbers </h1><br>""")
    for r in range(1, numofruns+1):
        Html_file.write("<table cellpadding=10><tr><td><h2>"+str(r)+": "+str(filenames1[r])+"</h2></td></tr>"+"<tr><td>Total sequences:</td><td>"+str(TSeqNums[r])+"</td></tr>"+"<tr><td>Passed sequences:</td><td>"+str(LSeqNums[r])+"</td></tr>"+"</h2></td></tr>"+"<tr><td>Barcode:</td><td>"+str(Barcodevalues[r])+"""</td></tr><tr><td style="background-color:lightgrey;"></td>""")
        for k in range(mink, maxk+1):
            Html_file.write("""<td style="background-color:lightgrey;"><b>K"""+str(k)+"</b></td>")
        Html_file.write("<tr><td>Total kmers</td>")
        for k in range(mink, maxk+1):
            Html_file.write("<td>"+str(numofkmers[r][k])+"</td>")
        Html_file.write("</tr><tr><td>Unique kmers</td>")
        for k in range(mink, maxk+1):
            Html_file.write("<td>"+str(numofuniquekmers[r][k])+" ("+str( round(((numofuniquekmers[r][k]/(4**k))*100),2)) +"%)"+"</td>")
        Html_file.write("</tr><tr><td>TFBS</td>")
        for k in range(mink, maxk+1):
            Html_file.write("<td>"+str(numoftfbs[r][k])+"</td>")
        Html_file.write("</tr><tr><td>Sequences w/ TFBS</td>")
        for k in range(mink, maxk+1):
            Html_file.write("<td>"+str(numoftfbsseq[r][k])+" ("+str(round(((numoftfbsseq[r][k]/LSeqNums[r])*100), 2))+"%)"+"</td>")
        Html_file.write("</tr><tr><td>Sequence bias</td>")
        for k in range(mink, maxk+1):
            Html_file.write("<td>"+str(round((sum(seqbias[r][k])/len(seqbias[r][k])),7))+"</td>")
        Html_file.write("</table>")

tablemaker()

Html_file.close()

print("DONE!")
