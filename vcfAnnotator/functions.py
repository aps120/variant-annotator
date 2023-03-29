import sys
import json
import time
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

"""functions for use with variant annotator"""

def find_repeat(s):
    '''find a repeating substring'''
    i = (s+s).find(s, 1, -1)
    return None if i == -1 else s[:i]
def dup_check(chrom1, pos1, ref1, alt1):
    '''check for a duplication'''
    dup = find_repeat(alt1[1:])
    if dup != None:
        hgvs = str(chrom1)+":g."+str(pos1+1)+"_"+str(pos1+len(dup))+"dup"+dup
    else:
        hgvs = None
    return hgvs

def vcf2hgvs(chrom1, pos1, ref1, alt1):
    '''convert each line of vcf to hgvs'''
    hgvs_start = str(chrom1) + ":g." + str(pos1)
    reflen = len(ref1)
    altlen = len(alt1)
    ## if it's a substitution
    if reflen == 1 and altlen == 1:
        hgvs = hgvs_start + ref1 + ">" + alt1
        var = "substitution"
    ## if it's an insertion
    elif reflen < altlen:
        ## first check for duplication
        hgvs = dup_check(chrom1, pos1, ref1, alt1)
        if hgvs != None:
            hgvs = hgvs
            var = "duplication"
        else:
            final_pos = str(pos1 + reflen)
            hgvs = hgvs_start + "_" + final_pos + "ins" + alt1
            var = "insertion"
    ## if it's a deletion
    elif reflen > altlen:
        final_pos = str(pos1 + reflen)
        hgvs = hgvs_start + "_" + final_pos + "del" + alt1
        var = "deletion"
    ##if it's an indel
    elif reflen == altlen:
        final_pos = str(pos1 + reflen)
        hgvs = hgvs_start + "_" + final_pos + "del" + ref1 + "ins" + alt1
        var = "indel"
    else:
        print("error converting variant at position %s to hgvs:" % pos1)
    return hgvs, var


class EnsemblRestClient(object):
    '''Class to work with Ensembl Rest API'''
    '''based on https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client'''
    ## it appeared that the reference genome was grch37 so I directed to there
    def __init__(self, server='http://grch37.rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_variants(self, hgvs):
        variants = self.perform_rest_action("/vep/human/hgvs/{0}".format(hgvs))
        if variants:
            return variants
        return None



def run(hgvs):
    '''grabs variants from Ensembl'''
    client = EnsemblRestClient()
    variant = client.get_variants(hgvs)
    return variant

def getVarInfo(hgvs):
    '''grabs variant info from the Ensembl Rest API'''
    r = run(hgvs)
    decoded = r[0]
    mydict = decoded
    ## get the effect
    effect = mydict['most_severe_consequence']
    ## get gene of variant
    key = "transcript_consequences"
    if mydict.get(key) == None:
        # if it's an intergenic region
        gene_symbol = 'intergenic'
    else:
        transcript = mydict['transcript_consequences'][0]
        gene_symbol = transcript["gene_symbol"]
    return effect, gene_symbol

