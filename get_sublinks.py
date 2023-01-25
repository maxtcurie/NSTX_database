from bs4 import BeautifulSoup
from urllib.request import Request, urlopen
import re

link="https://nstx.pppl.gov/nstx/SpecFit/"
def get_sublinks(link):
    req = Request("https://nstx.pppl.gov/nstx/SpecFit/")
    html_page = urlopen(req)

    soup = BeautifulSoup(html_page, "lxml")

    links = []
    for link in soup.findAll('a'):
        links.append(link.get('href'))
    return links


sublinks=get_sublinks(link)
#for the order by year

print(links)
