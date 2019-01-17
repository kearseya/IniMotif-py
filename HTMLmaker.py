from KmerKounter import identifier
from KmerKounter import numofruns
from KmerKounter import mink
from KmerKounter import maxk

def setidentity():
    global htmlname
    htmlname = identifier+'.html'
setidentity()


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

/* Style the header */
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
print(runheaders)


def results():
    html_strs = []
    for z in range(0, numofruns+1):
        html_strs.append({})
        for i in range(mink, maxk+1):
            html_strs[z].update({i:[]})

    for x in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            string = """<tr><td><p>K = """+str(k)+"""</p><img src="figures/logo_"""+str(identifier)+"_"+str(x)+"_"+str(k)+""".png"><img src="figures/hamdist_"""+str(identifier)+"_"+str(x)+"_"+str(k)+""".png"><img src="figures/pos_"""+str(identifier)+"_"+str(x)+"_"+str(k)+""".png"></td></tr>"""
            html_strs[x][k].append(string)

    return html_strs

html_strs = results()

def kmerfrequency():
    html_kmerfreq = []



def formatter():
    for r in range(1,numofruns+1):
        Html_file.write(runheaders[r])
        for k in range(mink, maxk+1):
            Html_file.write(html_strs[r][k][0])




formatter()


Html_file.close()

print("DONE!")
