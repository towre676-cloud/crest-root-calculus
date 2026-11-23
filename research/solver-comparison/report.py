"""Rebuild tables and paired CPU comparisons from the complete raw receipts."""
from collections import defaultdict
import json
from math import exp,log
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent


def analyze(study,case_ids=None):
    receipt = json.loads((HERE/(study+'-results.json')).read_text())
    protocol = json.loads((HERE/(study+'-results-protocol.json')).read_text())
    panel = json.loads((HERE/(study+'-panel.json')).read_text())
    if case_ids is not None:
        chosen = set(case_ids)
        panel = {**panel,'cases':[case for case in panel['cases'] if case['id'] in chosen]}
    chosen = {case['id'] for case in panel['cases']}
    if not receipt['completed']:
        raise ValueError('Wait for the complete frozen study.')
    groups = defaultdict(list)
    for row in receipt['runs']:
        if row['case'] in chosen:
            groups[row['case'],row['method']].append(row)
    medians = {}
    wall_medians = {}
    completions = {}
    for method in protocol['methods']:
        complete_cases = []
        for case in panel['cases']:
            rows = groups[case['id'],method]
            if len(rows)==2 and all(row['status']=='complete' for row in rows):
                medians[case['id'],method] = statistics.median(row['solve_cpu_seconds'] for row in rows)
                wall_medians[case['id'],method] = statistics.median(row['process_wall_seconds'] for row in rows)
                complete_cases.append(case['id'])
        completions[method] = {'completed_twice':len(complete_cases),'cases':len(panel['cases']),
                               'timeouts':sum(row['status']=='timeout' for row in receipt['runs'] if row['method']==method and row['case'] in chosen),
                               'errors':sum(row['status']=='error' for row in receipt['runs'] if row['method']==method and row['case'] in chosen),
                               'wrong_answers':sum(row['status']=='wrong_answer' for row in receipt['runs'] if row['method']==method and row['case'] in chosen)}
    reference = protocol['methods'][0]
    paired = {}
    for method in protocol['methods'][1:]:
        ratios = [medians[case['id'],method]/medians[case['id'],reference] for case in panel['cases']
                  if (case['id'],method) in medians and (case['id'],reference) in medians]
        wall_ratios = [wall_medians[case['id'],method]/wall_medians[case['id'],reference] for case in panel['cases']
                       if (case['id'],method) in medians and (case['id'],reference) in medians]
        if ratios:
            paired[method] = {'jointly_completed_cases':len(ratios),
                              'median_comparator_over_crest_cpu':statistics.median(ratios),
                              'median_crest_over_comparator_cpu':statistics.median(1/value for value in ratios),
                              'geometric_mean_comparator_over_crest_cpu':exp(statistics.mean(map(log,ratios))),
                              'median_comparator_over_crest_process_wall':statistics.median(wall_ratios),
                              'median_crest_over_comparator_process_wall':statistics.median(1/value for value in wall_ratios),
                              'geometric_mean_comparator_over_crest_process_wall':exp(statistics.mean(map(log,wall_ratios))),
                              'minimum_ratio':min(ratios),'maximum_ratio':max(ratios),
                              'crest_faster_cases':sum(ratio>1 for ratio in ratios)}
    return receipt,protocol,panel,groups,medians,{'completion':completions,'paired':paired}


def timings(rows):
    if len(rows)==2 and all(row['status']=='complete' for row in rows):
        return f"{statistics.median(row['solve_cpu_seconds'] for row in rows):.4f}"
    return ' / '.join('T' if row['status']=='timeout' else 'E' if row['status']=='error'
                      else f"W ({row['solve_cpu_seconds']:.4f})" if row['status']=='wrong_answer'
                      else f"{row['solve_cpu_seconds']:.4f}" for row in sorted(rows,key=lambda row:row['replicate']))


def main():
    sections = ["The three frozen panels contain 51 inputs and 312 planned runs in the published studies. Returned answers are checked exactly, and the map study uses a frozen verified reference. Wrong answers and timeouts remain in the recorded schedule. CPU seconds include each worker and its reaped native children; the adapter descriptions identify the work included.", "",
                "The quartic study resumed after 55 recorded runs without replacing any earlier measurement. Its input and adapter were frozen before the original start; imported research-module hashes were captured at resumption and verified at completion. The initial protocol stored the available affinity set, but not the selected worker CPU; the resumption records CPU eight. The later studies froze their mathematical dependencies before timing. The original 12-input counting panel and its 48 runs are preserved separately; the final counting panel retains those inputs and adds a predeclared three-input larger-field cohort, with every case retimed. A first map attempt stopped on an exact disagreement; its 11 recorded rows and original protocol are preserved separately. The published map study reran the unchanged input panel with verified references and a separately labelled corrected dense comparator. See the [protocol and reproduction notes](README.md).", ""]
    stats = {}
    for study,title in [('quartic','Marked quartic pencils'),('map','Constant-vertical hyperelliptic maps'),
                        ('frobenius','Hyperelliptic local zeta numerators')]:
        receipt,protocol,panel,groups,medians,summary = analyze(study)
        if study=='frobenius':
            summary['cohorts'] = {name:analyze(study,case_ids)[-1] for name,case_ids in panel['cohorts'].items()}
        stats[study] = summary
        sections += ['## '+title,'',f"Each input was run twice per method with a {protocol['per_process_wall_limit_seconds']:g}-second process limit. Entries below give median solver CPU seconds when both runs completed correctly. T records a timeout. W marks an incorrect answer and retains its CPU time in parentheses. E would mark a process error. A mixed entry preserves both outcomes in repetition order.",'']
        methods = protocol['methods']
        heading = '| Input | Exact outcome | '+' | '.join(methods)+' |'
        sections += [heading,'|---|---|'+'---:|'*len(methods)]
        for case in panel['cases']:
            name = case['id']
            if study=='quartic':
                count = receipt['supports'][name]['geometric_marks']
                outcome = f'{count} geometric pencil'+('s' if count!=1 else '')
            elif study=='map':
                output = receipt['outputs'][name]
                outcome = f"{output['geometric_maps']} geometric, {len(output['rational_maps'])} rational"
            else:
                outcome = 'degree '+str(2*case['genus'])+'; max field '+str(case['prime']**case['genus'])
            sections.append('| '+name+' | '+outcome+' | '+' | '.join(timings(groups[name,method]) for method in methods)+' |')
        sections += ['','| Method | Correct inputs completed twice | Timeout runs | Wrong-answer runs | Process-error runs |','|---|---:|---:|---:|---:|']
        for method,record in summary['completion'].items():
            sections.append(f"| {method} | {record['completed_twice']}/{record['cases']} | {record['timeouts']} | {record['wrong_answers']} | {record['errors']} |")
        sections += ['','Paired ratios use only inputs that both methods completed correctly in both repetitions. Each ratio divides the comparator’s median time by the project method’s median time on the same input. CPU ratios use solver CPU; process-wall ratios also include interpreter startup. Both use the same jointly completed inputs. A ratio above one means the comparator took more time. Timeouts and incorrect outputs contribute to their respective table columns and remain outside these finite ratios.','',
                     '| Comparator | Paired inputs | Median CPU ratio | Geometric mean CPU ratio | Median process-wall ratio | Project faster on CPU |',
                     '|---|---:|---:|---:|---:|---:|']
        for method,record in summary['paired'].items():
            sections.append(f"| {method} | {record['jointly_completed_cases']} | {record['median_comparator_over_crest_cpu']:.3f} | {record['geometric_mean_comparator_over_crest_cpu']:.3f} | {record['median_comparator_over_crest_process_wall']:.3f} | {record['crest_faster_cases']}/{record['jointly_completed_cases']} |")
        sections += ['']
        if study=='frobenius':
            sections += ['### Predeclared counting cohorts','',
                         'The small-prime cohort contains the twelve original inputs. The larger-field cohort contains three seeded inputs with maximum extension orders 1009² = 1,018,081, 97³ = 912,673 and 31⁴ = 923,521. The two-million-element table budget applies to the direct counter. Both cohorts use the same ten-second process and four-GiB memory limits. Cohort membership was fixed before the final timing study.','',
                         '| Cohort | Direct correct twice | PARI correct twice | Paired inputs | Median PARI/direct CPU ratio | Median PARI/direct process-wall ratio |',
                         '|---|---:|---:|---:|---:|---:|']
            for name,cohort in summary['cohorts'].items():
                first,second = cohort['completion']['crest-direct'],cohort['completion']['pari']
                pair = cohort['paired'].get('pari')
                values = (str(pair['jointly_completed_cases']),f"{pair['median_comparator_over_crest_cpu']:.4f}",f"{pair['median_comparator_over_crest_process_wall']:.4f}") if pair else ('0','—','—')
                sections.append(f"| {name} | {first['completed_twice']}/{first['cases']} | {second['completed_twice']}/{second['cases']} | "+' | '.join(values)+' |')
            sections += ['']
    (HERE/'results.md').write_text('# Measured solver comparisons\n\n'+'\n'.join(sections))
    (HERE/'summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(json.dumps(stats,indent=2))


if __name__=='__main__':
    main()
