"""Optional real-browser inspection of the offline workshop (Playwright and Chromium)."""
import argparse
import json
from pathlib import Path
import sys
import threading
from http.server import ThreadingHTTPServer
from playwright.sync_api import sync_playwright
from workshop import CATALOGUE, Handler

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chromium',type=Path)
    parser.add_argument('--exact',action='store_true',help='also run the optional Part XXX engine')
    args=parser.parse_args()
    with ThreadingHTTPServer(('127.0.0.1',0),Handler) as server,sync_playwright() as play:
        threading.Thread(target=server.serve_forever,daemon=True).start()
        url=f'http://127.0.0.1:{server.server_address[1]}/'
        browser=play.chromium.launch(executable_path=str(args.chromium) if args.chromium else None,
                                     args=['--no-sandbox'])
        page=browser.new_page(viewport={'width':1200,'height':850})
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(url);page.wait_for_function('window.location && document.querySelector("#answer")')
        for key,info in CATALOGUE.items():
            page.locator(f'#apps a[href="#{key}"]').click()
            page.wait_for_function('(s)=>document.querySelector("#domain").textContent===s',arg=info['name'])
            for i,case in enumerate(info['cases']):
                if case['file'].startswith('original-') and not args.exact: continue
                if i: page.locator('#cases').select_option(case['file'])
                page.wait_for_function('document.querySelector("#activity").textContent==="Calculation complete."',
                                      timeout=80000)
                assert page.locator('#out pre').count()==1
                result=json.loads(page.locator('#out pre').text_content())
                assert result
                print(f"{key}: {case['file']} OK",flush=True)
            guide=page.request.get(url+'guide/'+key)
            assert guide.ok and '<main>' in guide.text()
        page.locator('#apps a[href="#wafer"]').click()
        page.wait_for_selector('.follow input')
        page.locator('.follow input').fill('0.27735')
        page.get_by_role('button',name='Apply measurement').click()
        page.wait_for_function('document.querySelector("#out pre")?.textContent.includes("within_limit")')
        with page.expect_download() as event:page.get_by_role('button',name='Save calculation record').click()
        record=json.loads(Path(event.value.path()).read_text())
        assert record['result']['surviving_profiles']==3 and '32' in record['input']['observations']
        page.set_viewport_size({'width':360,'height':780})
        page.wait_for_function('document.documentElement.scrollWidth<=window.innerWidth')
        assert page.locator('#apps a').count()==7
        assert not errors,errors
        browser.close();server.shutdown()
        print('All requested browser cases, measured follow-up, export and mobile layout passed.')

if __name__=='__main__':main()
