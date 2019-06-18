
import requests

url = "http://rest.kegg.jp/get/hsa00010/kgml"

import urllib.request as urllib2

from xml.dom import minidom

dom = minidom.parse(urllib2.urlopen(url))
