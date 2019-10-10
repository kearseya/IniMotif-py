import requests
from xml.etree.ElementTree import fromstring, ElementTree
import re
import urllib.request

def getXML(search, dataType):
    # download an xml file for the specified search terms	

    # build the url for the query and download the xml file
    build_url = {"query": search,
                 "result": dataType,
                 "length": "1",
                 "download": "xml",
                 }
    
    response = requests.get("https://www.ebi.ac.uk/ena/data/search", params=build_url)
    xmlout = response.content

    return xmlout


def parseXMLgetFTP(xmlout):
    # parse the xml for http links which contain information on the fastq files
    # open the http links and write the result to file
    
    # create element tree object
    tree = ElementTree(fromstring(xmlout))
    # get root element
    root = tree.getroot()
    
    httplinks = []
    # iterate xml file for http links
    for item in root.iter("ID"):
        if item.text.startswith("http://") and item.text.endswith("fastq_bytes"):
            httplinks.append(item.text)

    ftpinfo = []
    for url in httplinks:
        response = requests.get(url)
        ftpinfo.append(response.content)

    return ftpinfo


def parseFTPgetFASTQ(ftpinfo):
    # parse for the ftp links and download
    
    # convert ftpinfo (bytes) to string
    ftpstring = [x.decode() for x in ftpinfo]
    # use regex to search ftp links
    regexFTP = re.compile("ftp.sra(.*).fastq.gz")
    
    for item in ftpstring:
        match = regexFTP.search(item)
        result = match.group(0)
        for elem in result.split(";", 2):
            filename = elem[elem.rfind("/")+1:]
            ftplink = "ftp://" + elem
            urllib.request.urlretrieve(ftplink, filename)


xmlout = getXML("PRJNA385825", "READ_STUDY")
ftpinfo = parseXMLgetFTP(xmlout)
parseFTPgetFASTQ(ftpinfo)
