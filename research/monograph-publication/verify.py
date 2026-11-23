"""Check the archived monograph or retrieve and check its public edition."""
import argparse
from collections import Counter
from hashlib import sha256
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
URL="https://towre676-cloud.github.io/crest-root-calculus/"
SHA256="978a72fe6d31b629aa18ca72cd4aa3649fe8a6191e5d8656693861aa047961a0"


class Contents(HTMLParser):
    def __init__(self):
        super().__init__()
        self.sections=[]
        self.fragments=set()
        self.anchors=Counter()
        self.heading=None

    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if "id" in attrs:
            self.anchors[attrs["id"]]+=1
        if tag=="a" and attrs.get("href","").startswith("#"):
            self.fragments.add(attrs["href"][1:])
        if tag=="h2":
            self.heading=[]

    def handle_data(self,text):
        if self.heading is not None:
            self.heading.append(text)

    def handle_endtag(self,tag):
        if tag=="h2" and self.heading is not None:
            title="".join(self.heading)
            number=re.match(r"^(\d+)\.\s",title)
            if number:
                self.sections.append((int(number[1]),title))
            self.heading=None


def verify(data):
    digest=sha256(data).hexdigest()
    if digest!=SHA256:
        raise ValueError("The edition differs from crest/index.html in the archived Part XXX release")
    contents=Contents()
    contents.feed(data.decode("utf-8"))
    if [n for n,t in contents.sections]!=list(range(1,149)):
        raise ValueError("Unexpected numbered-section sequence")
    if max(contents.anchors.values())!=1 or contents.fragments-set(contents.anchors):
        raise ValueError("Broken internal anchor")
    return {"bytes":len(data),"sha256":digest,"numbered_sections":len(contents.sections),
            "first":contents.sections[0][1],"last":contents.sections[-1][1]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live",action="store_true",help="compare the HTTPS response with the archived bytes")
    args=parser.parse_args()
    result=verify((HERE/"index.html").read_bytes())
    if args.live:
        with urlopen(Request(URL,headers={"User-Agent":"Crest-Root edition verifier"}),timeout=30) as response:
            if response.status!=200 or response.headers.get_content_type()!="text/html" or response.geturl()!=URL:
                raise ValueError("Unexpected status, media type or public address")
            if verify(response.read())!=result:
                raise ValueError("Published edition differs from local edition")
    print(json.dumps(result,indent=2))
