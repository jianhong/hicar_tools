#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
import os

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    "nf-core/hicar": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
}
results = OrderedDict()
results["nf-core/hicar"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Search software.version.txt
version_files = [x for x in os.listdir('.') if x.endswith('.version.txt')]
for version_file in version_files:

    software = version_file.replace('.version.txt','')

    with open(version_file) as fin:
        version = fin.read().strip()

    results[software] = version

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/hicar Software Versions'
section_href: 'https://github.com/nf-core/hicar'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
