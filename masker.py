import re
#from Bio import SeqIO
import os
import shutil
import fileinput

from itertools import chain, combinations, product

inputlist = []

firsttime = True

try:
    from tkinter import *

    variableentrylist = []
    variablevaluelist = []

    currentshown = 1
    currenthidden = 0

    maxmasks = 1

    def changeform(event):
        global variableentrylist
        global variablevaluelist
        global currentshown
        global firsttime
        global maxmasks
        masktypelist = ["Repeat", "Motif"]
        nummasks = int(numofmasksentry.get())

        if nummasks*4 > len(variableentrylist):
            if firsttime == True:
                maxmasks = 1

            for x in range((maxmasks+1), (nummasks+1)):
                variableentrylist.append("typeent"+str(x))
                variablevaluelist.append("typeval"+str(x))
                variablevaluelist[(x*4)-4] = StringVar()
                variableentrylist[(x*4)-4] = OptionMenu(maskinputframe, variablevaluelist[(x*4)-4], "Repeat", "Motif")
                variablevaluelist[(x*4)-4].set(masktypelist[0])
                variableentrylist[(x*4)-4].grid(row=0, column=x)
                variableentrylist.append("unitent"+str(x))
                variablevaluelist.append("unitval"+str(x))
                variablevaluelist[(x*4)-3] = StringVar()
                variableentrylist[(x*4)-3] = Entry(maskinputframe, textvariable=variablevaluelist[(x*4)-3], width=10)
                variableentrylist[(x*4)-3].grid(row=1, column=x)
                variableentrylist.append("revent"+str(x))
                variablevaluelist.append("revval"+str(x))
                variablevaluelist[(x*4)-2] = BooleanVar()
                variableentrylist[(x*4)-2] = Checkbutton(maskinputframe, variable=variablevaluelist[(x*4)-2])
                variableentrylist[(x*4)-2].grid(row=2, column=x)
                variableentrylist.append("rephament"+str(x))
                variablevaluelist.append("rephamval"+str(x))
                variablevaluelist[(x*4)-1] = IntVar()
                variableentrylist[(x*4)-1] = Entry(maskinputframe, textvariable=variablevaluelist[(x*4)-1], width=10)
                variableentrylist[(x*4)-1].delete(0)
                variableentrylist[(x*4)-1].grid(row=3, column=x)
                if x == nummasks:
                    firsttime = False
                    maxmasks = nummasks
                    currentshown = nummasks

        if nummasks > currentshown:
            for j in range(currentshown+1, nummasks+2):
                variableentrylist[((j-1)*4)-4].grid(row=0, column=j)
                variableentrylist[((j-1)*4)-3].grid(row=1, column=j)
                variableentrylist[((j-1)*4)-2].grid(row=2, column=j)
                variableentrylist[((j-1)*4)-1].grid(row=3, column=j)
                if j == nummasks+1:
                    currentshown = nummasks

        if nummasks < currentshown:
            for i in range(nummasks+1, int(len(variableentrylist)/4)+1):
                variableentrylist[(i*4)-4].grid_forget()
                variableentrylist[(i*4)-3].grid_forget()
                variableentrylist[(i*4)-2].grid_forget()
                variableentrylist[(i*4)-1].grid_forget()
                if i == int(len(variableentrylist)/4):
                    currentshown = nummasks

    def makeinput():
        global file
        file = str(filelocationentry.get())
        global numberofmasks
        numberofmasks = int(numofmasksentry.get())
        global inputlist
        for x in range(1, numberofmasks+1):
            inputlist.append(str(variablevaluelist[(x*4)-4].get()))
            inputlist.append(str(variablevaluelist[(x*4)-3].get()))
            rev = variablevaluelist[(x*4)-2].get()
            if rev == True:
                inputlist.append("yes")
            if rev == False:
                inputlist.append("no")
            inputlist.append(int(variablevaluelist[(x*4)-1].get()))


    window = Tk()

    titleframe = Frame(window)
    titleframe.pack(side=TOP)

    title = Label(titleframe, text="Masker")
    title.pack(side=TOP)

    submitframe = Frame(window)
    submitframe.pack(side=BOTTOM, anchor="n", pady=10)

    maskinputframe = Frame(window)
    maskinputframe.pack(side=BOTTOM, anchor="n", padx=20, pady=20)

    inputframe = Frame(window)
    inputframe.pack(side=LEFT)

    pictureframe = Frame(window)
    pictureframe.pack(side=RIGHT)
    maskerimage = PhotoImage(file='figures/GUIgraphics/masker.png')
    maskerimagelabel = Label(pictureframe, image=maskerimage)
    maskerimagelabel.grid(padx=20, pady=20)

    filelocationlabel = Label(inputframe, text="File with path: ")
    filelocationlabel.grid(row=0, column=0, padx=10, sticky="e")
    filelocationentry = Entry(inputframe, textvariable=StringVar())
    filelocationentry.grid(row=0, column=1)

    numofmaskslabel = Label(inputframe, text="Number of masks: ")
    numofmaskslabel.grid(row=1, column=0, padx=10, sticky="e")
    numofmasksentry = Entry(inputframe, textvariable=IntVar())
    numofmasksentry.grid(row=1, column=1)
    numofmasksentry.delete(0)
    numofmasksentry.insert(0, int(1))
    numofmasksentry.bind("<FocusOut>", changeform)
    numofmasksentry.bind("<Return>", changeform)

    typelabel = Label(maskinputframe, text="Type of masks: ")
    typelabel.grid(row=0, column=0, padx=10, sticky="e")
    unitlabel = Label(maskinputframe, text="Unit string: ")
    unitlabel.grid(row=1, column=0, padx=10, sticky="e")
    revwantlabel = Label(maskinputframe, text="Reverse compliment: ")
    revwantlabel.grid(row=2, column=0, padx=10, sticky="e")
    minrephamlabel = Label(maskinputframe, text="Min repeats/Mutations: ")
    minrephamlabel.grid(row=3, column=0, padx=10, sticky="e")

    variableentrylist.append("typeent1")
    variablevaluelist.append("typeval1")
    variablevaluelist[0] = StringVar()
    variableentrylist[0] = OptionMenu(maskinputframe, variablevaluelist[0], "Repeat", "Motif")
    variablevaluelist[0].set("Repeat")
    variableentrylist[0].grid(row=0, column=1)
    variableentrylist.append("unitent1")
    variablevaluelist.append("unitval1")
    variablevaluelist[1] = StringVar()
    variableentrylist[1] = Entry(maskinputframe, textvariable=variablevaluelist[1], width=10)
    variableentrylist[1].grid(row=1, column=1)
    variableentrylist.append("revent1")
    variablevaluelist.append("revval1")
    variablevaluelist[2] = BooleanVar()
    variableentrylist[2] = Checkbutton(maskinputframe, variable=variablevaluelist[2])
    variableentrylist[2].grid(row=2, column=1)
    variableentrylist.append("rephament1")
    variablevaluelist.append("rephamval1")
    variablevaluelist[3] = IntVar()
    variableentrylist[3] = Entry(maskinputframe, textvariable=variablevaluelist[3], width=10)
    variableentrylist[3].delete(0)
    variableentrylist[3].grid(row=3, column=1)

    submitbutton = Button(submitframe, text="SUBMIT", padx=1, pady=1, command=makeinput)
    submitbutton.pack(anchor="s")

    window.mainloop()

except:
    def clinputbuffer():
        global file
        global numberofmasks
        global inputlist
        inputlist =[]
        file = str(input("File to be masked (with path): "))
        numberofmasks = int(input("Number of masks: "))
        for x in range(1, numberofmasks+1):
            type = str(input("Mask type (repeat/motif): "))
            inputlist.append(type)
            if type in ["repeat", "r", "rep"]:
                unit = str(input("Unit string: "))
                inputlist.append(unit)
                revwant = str(input("Mask revcomp of unit?: "))
                inputlist.append(revwant)
                min_rep = int(input("Minimum repeats: "))
                inputlist.append(min_rep)
            if type in ["motif", "m", "mot"]:
                unit = str(input("Unit string: "))
                inputlist.append(unit)
                revwant = str(input("Mask revcomp of unit?: "))
                inputlist.append(revwant)
                allowham = int(input("Num mutations: "))
                inputlist.append(allowham)

    clinputbuffer()



revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev

def hamming_circle(s, n, alphabet):
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)

def hamming_ball(s, n, alphabet):
    allowedkmers = set(chain.from_iterable(hamming_circle(s, i, alphabet) for i in range(n + 1)))
    return allowedkmers

def myrepl(m):
    return ('N'*len(m.group(0)))

#pat = re.compile(f'({unit_str}){{{n_min_unit},{n_max_unit}}}')
#out_str = re.sub(pat, myrepl, instr)
max_rep = ''


def masker():
    outfilename = str(os.path.dirname(file))+"/masked_"+str(os.path.basename(file))
    shutil.copy(file, str(os.path.dirname(file))+"/masked_"+str(os.path.basename(file)))
    for n in range(1, numberofmasks+1):
        typeofmask = inputlist[(n*4)-4]
        if typeofmask in ["Repeat", "repeat", "r", "R", "rep"]:
            unit = inputlist[(n*4)-3]
            revwant = inputlist[(n*4)-2]
            min_rep = inputlist[(n*4)-1]
            funit = {unit}
            pat = re.compile(f'({unit}){{{min_rep},{max_rep}}}')
            frunit = set()
            if revwant in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
                for i in funit:
                    frunit.add(i)
                    frunit.add(revComp(i))
                for unit in frunit:
                    pat = re.compile(f'({unit}){{{min_rep},{max_rep}}}')
                    for line in fileinput.input([outfilename], inplace=True):
                        print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))
            else:
                pat = re.compile(f'({unit}){{{min_rep},{max_rep}}}')
                for line in fileinput.input([outfilename], inplace=True):
                    print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))

        if typeofmask in ["Motif", "motif", "m", "M", "mot"]:
            unit = inputlist[(n*4)-3]
            min_rep = 1
            revwant = inputlist[(n*4)-2]
            allowham = inputlist[(n*4)-1]
            maskmotifs = hamming_ball(unit, allowham, "ATGC")
            rmaskmotifs = set()
            if revwant in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
                for i in maskmotifs:
                    rmaskmotifs.add(i)
                    rmaskmotifs.add(revComp(i))
                for kmer in rmaskmotifs:
                    pat = re.compile(f'({kmer}){{{min_rep},{max_rep}}}')
                    for line in fileinput.input([outfilename], inplace=True):
                        print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))
            else:
                for kmer in maskmotifs:
                    pat = re.compile(f'({kmer}){{{min_rep},{max_rep}}}')
                    for line in fileinput.input([outfilename], inplace=True):
                        print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))

masker()
