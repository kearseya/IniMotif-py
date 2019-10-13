from KmerKounter import identifier
from KmerKounter import numofruns
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import filenames1
from KmerKounter import startround
from KmerKounter import multiround
if multiround == True:
    from KmerKounter import files

from KmerKounter import totaldict

from WebLogoMod import nmotifs
from WebLogoMod import countdict

from PositionBias import TSeqNums
from PositionBias import LSeqNums
#from PositionBias import Barcodevalues
from KmerKounter import barcodeprimers53
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
  background-image: linear-gradient(to bottom, #051937, #004d7a, #008793, #00bf72, #a8eb12);
}

p {
  font-size: 100px;
}

h1 {
  font-size: 200px;
  background-color: none;
  color: white;
}

h2 {
  font-size: 35px;
}

th {
  color: lightgrey;
  background-color: none;
}

td {
  background-color: white;
}

header {
  background-color: none;
  padding: 10px;
  text-align: center;
  font-size: 35px;
  color: white;

}
</style>
</head>
<header>
<h1 style="font-size: 600px;"> IniMotif </h1>
<h2 style="font-size: 100px">A pipeline for motif discovery!</h2>
</header>
<br>

"""

logoordernames = ["", "First", "Second", "Third", "Fourth", "Fifth", "Sixth"]

Html_file.write(title)

def runheader():
    runheaders = []
    for x in range(0, numofruns+1):
        string= """<h1 align = "middle" colour = "white";>Run """+str(x+(startround-1))+"</h1>"
        runheaders.append(string)
    return runheaders

runheaders = runheader()

# style="background-color: lightgrey;
def results():
    html_strs = []
    for z in range(0, numofruns+1):
        html_strs.append({})
        for i in range(mink, maxk+1):
            html_strs[z][i] ={}
            for n in range(1, nmotifs+1):
                html_strs[z][i][n] = []

    for x in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            #logos = "<tr>"
            #kval = """<p style="color: white; padding: 30px;" align = "middle"><b>K = """+str(k)+"</b></p>"""
            for n in range(1, nmotifs+1):
                logos = """<tr><td colspan = 2 align = "middle"><img src="figures/"""+str(identifier)+"""/logos/logo_"""+str(identifier)+"_"+str(x+(startround-1))+"_"+str(k)+"_"+str(n)+""".png"></td></tr>"""
                if nmotifs == 1:
                    logosusage = """<tr><td colspan = 2 align = "middle"><p>"""+"Motif uses "+str(round((countdict[x][k][n]/totaldict[x][k])*100, 2))+"% of kmers"+"</p></td></tr>"
                else:
                    logosusage = """<tr><td colspan = 2 align = "middle"><p>"""+str(logoordernames[n])+" motif uses "+str(round((countdict[x][k][n]/totaldict[x][k])*100, 2))+"% of kmers"+"</p></td></tr>"
                ham = """<tr><td><img src="figures/"""+str(identifier)+"""/hamming_distance/hamdist_"""+str(identifier)+"_"+str(x+(startround-1))+"_"+str(k)+"_"+str(n)+""".png" width=2000px hight=1200px></td>"""
                pos = """<td><img src="figures/"""+str(identifier)+"""/position_bias/pos_"""+str(identifier)+"_"+str(x+(startround-1))+"_"+str(k)+"_"+str(n)+""".png" width=2000px hight=1200px></td></tr>"""
                string = """<table align = "center";>"""+logos+logosusage+ham+pos+"</table>"
                html_strs[x][k][n].append(string)

    return html_strs

html_strs = results()

#else:
    #logos += """<td><img src="figures/"""+str(identifier)+"""/logos/logo_"""+str(identifier)+"_"+str(x+(startround-1))+"_"+str(k)+"_"+str(n)+""".png"></td>""
#"""<p style="color: white; padding: 30px;" width=3000px hight=1800px align = "middle">K = """+str(i)+"""</p><img src="figures/"""+str(identifier)+"""/kmer_frequency/kmerfreq_"""+str(identifier)+"_"+str(k)+""".png" width=3000px height=3000px>"""
def kmerfrequency():
    html_kmerfreq = {}
    twos = 0
    #pair = """<table align = "center";>"""
    for k in range(mink, maxk+1):
        html_kmerfreq.update({k:[]})
        if twos % 2 == 0:
            first = """<tr><td><img src="figures/"""+str(identifier)+"""/kmer_frequency/kmerfreq_"""+str(identifier)+"_"+str(k)+""".png" width=3000px height=3000px></td>"""
            html_kmerfreq[k].append(first)
        if twos % 2 == 1:
            second = """<td><img src="figures/"""+str(identifier)+"""/kmer_frequency/kmerfreq_"""+str(identifier)+"_"+str(k)+""".png" width=3000px height=3000px></td></tr>"""
            html_kmerfreq[k].append(second)
        #print(html_kmerfreq)
        #html_kmerfreq.update({k:[]})
        twos += 1
    return html_kmerfreq

html_kmerfreq = kmerfrequency()
#print(html_kmerfreq)


def formatter():
    for r in range(1,numofruns+1):
        Html_file.write(runheaders[r])
        for k in range(mink, maxk+1):
            Html_file.write("""<p style="color: white; padding: 30px;" align = "middle"><b>K = """+str(k)+"</b></p>")
            for n in range(1, nmotifs+1):
                Html_file.write(html_strs[r][k][n][0])

def kmerfreqformatter():
    Html_file.write("""<h1 style="background-color: #3c3c3c; padding: 0px; color: white;" align = "middle"> Kmer frequency </h1>""")
    Html_file.write("""<table align = "center";>""")
    twos = 0
    for k in range(mink, maxk+1):
        Html_file.write(html_kmerfreq[k][0])
        twos += 1
        if k == maxk:
            if twos % 2 == 1:
                Html_file.write("""<td><background color = "white" width=3000px height=3000px></td></tr>""")
    Html_file.write("</table>")


formatter()

kmerfreqformatter()

#for _ in range(mink, maxk+1):
    #Html_file.write("<td></td>")

def tablemaker():
    Html_file.write("""<h1 style="background-color: #3c3c3c; padding: 20px; color: white;" align = "middle"> Numbers </h1><br>""")
    for r in range(1, numofruns+1):
        if multiround == False:
            Html_file.write("<table cellpadding=10><tr><td colspan = "+str((maxk-mink)+2)+""" align = "middle"><h2>"""+str(r+(startround-1))+": "+str(filenames1[r])+"</h2></td>")
            Html_file.write("</tr><tr><td>Total sequences:</td><td colspan = "+str((maxk-mink)+1)+""" align = "left">"""+str(TSeqNums[r])+"</td>")
            Html_file.write("</tr><tr><td>Passed sequences:</td><td colspan = "+str((maxk-mink)+1)+""" align = "left">"""+str(LSeqNums[r])+"</td>")
            Html_file.write("</tr></h2></td></tr>"+"<tr><td>Barcode:</td><td colspan = "+str((maxk-mink)+1)+""" align = "left">5': """+str(barcodeprimers53[r][0])+" ("+str(len(barcodeprimers53[r][0]))+"),   3': "+str(barcodeprimers53[r][1])+" ("+str(len(barcodeprimers53[r][1]))+")</td>")
            Html_file.write("""</tr><tr><td style="background-color:lightgrey;"></td>""")
        if multiround == True:
            Html_file.write("<table cellpadding=10><tr><td colspan = "+str((maxk-mink)+2)+""" align = "middle"><h2>"""+str(r+(startround-1))+": "+str(files[0])+"</h2></td>")
            Html_file.write("</tr><tr><td>Total sequences:</td><td colspan = "+str((maxk-mink)+1)+""" align = "left">"""+str(TSeqNums[r])+"</td>")
            Html_file.write("</tr><tr><td>Passed sequences:</td><td colspan = "+str((maxk-mink)+1)+""" align = "left">"""+str(LSeqNums[r])+"</td>")
            Html_file.write("</tr></h2></td></tr>"+"<tr><td>Barcode:</td><td colspan = "+str((maxk-mink)+1)+""" align = "left">5': """+str(barcodeprimers53[r][0])+" ("+str(len(barcodeprimers53[r][0]))+"),   3': "+str(barcodeprimers53[r][1])+" ("+str(len(barcodeprimers53[r][1]))+")</td>")
            Html_file.write("""</tr><tr><td style="background-color:lightgrey;"></td>""")
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
