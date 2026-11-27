from pathlib import Path
from playwright.sync_api import sync_playwright
import argparse,json,tempfile
parser=argparse.ArgumentParser(description="Run the Linkwright cartoner browser acceptance scenarios.")
parser.add_argument("--chromium", help="Path to a Chromium executable")
parser.add_argument("--url", help="Deployed trial URL; defaults to the offline file")
args=parser.parse_args()
base=Path(__file__).resolve().parent
url=args.url or (base/'index.html').as_uri()
with sync_playwright() as pw, tempfile.TemporaryDirectory() as td:
 b=pw.chromium.launch(executable_path=args.chromium,args=['--no-sandbox'])
 pg=b.new_page(viewport={'width':1320,'height':920},accept_downloads=True)
 errors=[]
 pg.on('pageerror',lambda exc:errors.append(str(exc)))
 pg.goto(url,wait_until='networkidle',timeout=45000)
 assert "Measured current" in pg.locator('#job').inner_text()
 assert pg.locator('#axisLabel').text_content()=='Elapsed time'
 assert '21 mm' in pg.locator('#loops').inner_text()
 pg.locator('#case').select_option('opener')
 assert pg.locator('#job').inner_text()=='Carton opening and transfer'
 assert pg.locator('#overall').inner_text()=='More evidence needed'
 assert pg.locator('#graph svg circle').count()>=3
 pg.locator('#cycle').select_option(index=1);assert 'More evidence needed' not in pg.locator('#checks').inner_text()
 pg.locator('#channel').select_option(index=1); pg.locator('#view').select_option('all')
 pg.locator('#case').select_option('folder');assert pg.locator('#overall').inner_text()=='Requirement violated'
 pg.locator('#case').select_option('reach');assert '-0.00006' in pg.locator('#loops').inner_text()
 pg.locator('#inputs').click()
 conf=json.loads((base/'examples/opener.json').read_text());conf['name']='Historical bench run';conf['source']='machine';conf['setup']='bench-r5'
 f=Path(td)/'my.json';f.write_text(json.dumps(conf))
 csv=Path(td)/'my.csv';csv.write_text((base/'examples/opener.csv').read_text().replace('opener-r1','bench-r5'))
 pg.locator('#setupFile').set_input_files(f);pg.wait_for_function("document.querySelector('#setupText').value.includes('bench-r5')")
 pg.locator('#traceFile').set_input_files(csv);pg.wait_for_function("document.querySelector('#traceText').value.includes('bench-r5')")
 pg.locator('#run').click();assert pg.locator('#job').inner_text()=='Historical bench run'
 pg.get_by_text('Trial comparison',exact=True).click();pg.locator('#baseline').select_option('Passed');pg.locator('#notes').fill('Run QA-4')
 with pg.expect_download() as d:pg.locator('#export').click()
 file=Path(td)/'record.json';d.value.save_as(file); rec=json.loads(file.read_text());assert rec['comparison']['baseline']=='Passed'
 pg.locator('#case').select_option('folder')
 pg.locator('#recordFile').set_input_files(file);pg.wait_for_function("document.querySelector('#notice').textContent.includes('match the replay')")
 assert pg.locator('#notes').input_value()=='Run QA-4'
 rec['result']['status']='meets';file.write_text(json.dumps(rec));pg.locator('#recordFile').set_input_files(file)
 pg.wait_for_function("document.querySelector('#notice').textContent.includes('disagrees')")
 assert pg.locator('#report').is_hidden() and pg.locator('#export').is_disabled()
 pg.locator('#case').select_option('opener');pg.locator('#inputs').click();pg.locator('#traceText').fill('bad');pg.locator('#run').click()
 assert pg.locator('#report').is_hidden()
 pg.locator('#reset').click();assert pg.locator('#report').is_visible()
 pg.locator('#inputs').click()
 pg.set_viewport_size({'width':390,'height':844})
 assert pg.evaluate('document.documentElement.scrollWidth<=window.innerWidth')
 pg.locator('#case').select_option('folder');pg.locator('#cycle').select_option(index=1)
 pg.get_by_text('Trial guide',exact=True).click();assert pg.locator('h1').inner_text()=='Run a cartoner trial'
 pg.get_by_text('Open trial',exact=True).first.click()
 pg.get_by_text('Measured evidence',exact=True).click()
 assert pg.locator('#peak').inner_text()=='48.42 mA'
 assert pg.locator('#wave svg circle').count()==3
 assert pg.locator('#piade tr').count()==5 and pg.locator('#comopi tr').count()==8 and pg.locator('#alpi tr').count()==20
 assert '444,834' in pg.locator('#alarmSummary').inner_text()
 assert pg.locator('#group option').count()==20
 assert pg.evaluate('document.documentElement.scrollWidth<=window.innerWidth')
 pg.get_by_text('Held-out results and original files').click()
 assert '226,431' in pg.locator('#pooled').inner_text()
 pg.locator('#group').select_option('honey/Sine_0')
 assert pg.locator('#runs a').count()==3
 assert pg.locator('#runs a').first.get_attribute('href').startswith('https://darus.uni-stuttgart.de/')
 pg.get_by_text('Estimate the value of a trial on one machine').click()
 pg.locator('#share').fill('10');pg.locator('#value').fill('100');pg.locator('#cost').fill('100')
 assert 'net for this observation period' in pg.locator('#economics').inner_text()
 pg.locator('#share').fill('101');assert 'Enter a recoverable share' in pg.locator('#economics').inner_text()
 pg.locator('#sparseRead').click()
 assert pg.locator('#overall').inner_text()=='More evidence needed'
 pg.locator('#inputs').click();pg.get_by_text('Import a native four-bar servo recording').click()
 pg.locator('#nativeLimit').fill('10');pg.locator('#nativeError').fill('0')
 native=Path(td)/'example.txt'
 native.write_text('% Time [s], Goal Current [mA], Pres. Current [mA], Pres. Voltage [V], Pres. Position [rad], Pres. Velocity [rad/s]\n0.001,1,9,12.1,1,-3\n0.005,1,11,12.1,0,-4')
 pg.locator('#nativeFile').set_input_files(native)
 pg.wait_for_function("document.querySelector('#job').textContent==='example.txt'")
 assert pg.locator('#overall').inner_text()=='Requirement violated'
 pg.get_by_text('Trial comparison',exact=True).click()
 pg.locator('#baseline').select_option('Passed');pg.locator('#observed').select_option('Confirmed deviation')
 with pg.expect_download() as new_record:pg.locator('#export').click()
 paired=Path(td)/'paired.json';new_record.value.save_as(paired)
 pg.get_by_text('Measured evidence',exact=True).click();pg.locator('#auditFiles').set_input_files(paired)
 pg.wait_for_function("document.querySelector('#auditStatus').textContent.startsWith('1 paired')")
 assert pg.locator('#auditTable tr').first.locator('td').nth(2).inner_text()=='1'
 with pg.expect_download() as audit_download:pg.locator('#auditSave').click()
 audit_file=Path(td)/'audit.json';audit_download.value.save_as(audit_file)
 assert json.loads(audit_file.read_text())['counts']['trialOnly']==1
 assert not errors,errors
 print('PASS: stroke checks, measured holdouts, native import, sparse replay, production logs, economics, mobile, record export and paired audit')
