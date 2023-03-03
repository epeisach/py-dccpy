#===========================================================
# this module is to handle all the TLS related problems
#===========================================================



import os,sys,  shutil, math
from  config import *
from util import *

##########################################################
def check_tls(pdb, id):
    '''
    id = 0, pdb is a file; id=1, pdb is a list
    '''

    fp=pdb
    if id==0 and os.path.exists(pdb) : fp = open(pdb,'r').readlines()
    tls, pdbn = remove_pdb_item(fp, 1, 'tls')
    if len(tls)<15 : return

    chain, ch_range = chain_res_list(pdbn, 1)
    ntls,  prog = '?' , '?'
    all_range=[] #[['ch', ni, nf, ntls] ...] ni/nf=first/last residue number
    prog=os.popen('grep "REMARK   3   PROGRAM     :" -m 1  %s' %pdb).read().split(':')[1].strip()
    
    for i, x in enumerate(tls):
            
        if 'TLS GROUP :' in x [:25]:
            ntls=x.split(':')[1].strip()
            if not is_number(ntls):
                t1='Error: TLS group (=%s) is not a number \n' %(ntls)
                perror(t1)
                continue
            ntls=int(ntls)
            
        elif (('REMARK   3    RESIDUE RANGE :' in x[:29]) or
              ('REMARK   3    SELECTION:' in x[:29] and
               'buster' in prog.lower() and '{' not in x[20:])):
           
            if 'SELECTION:' in x:   #buster changed set to  SELECTION
                xx = 'REMARK   3    RESIDUE RANGE :' + x[24:48]
                t = parse_tls_range_refmac(xx, ntls)
            else:   #refmac 
                t = parse_tls_range_refmac(x, ntls)
            t.append(ntls)
            all_range.append(t)
#            print 'from refmac=',  t 

            
        elif 'REMARK   3    SELECTION:' in x[:29] and 'phenix' in prog.lower():  #phenix,
            
            tls_str = tls_selection_one_string(tls, i, 25)
            t = parse_tls_range_phenix(tls_str, ntls, chain, ch_range,fp)
           
            for y in t:
                y.append(ntls)
                all_range.append(y)
            #print 'from phenix=', t
         
        elif ('REMARK   3    SET :' in x [:20]  #buster in original format
              or 'SELECTION: {' in x[11:29]) :  #buster in modified format
            tls_str = tls_selection_one_string(tls, i, 20)
            t = parse_tls_range_buster(tls_str, ntls, chain, ch_range)
            for y in t:
                y.append(ntls)
                all_range.append(y)
           # print 'from buster=', t
            
#    print all_range
    if not all_range : return
    for x in all_range:
        if len(x)<4:  continue
        ch,n1,n2,nt=x[0], x[1], x[2], x[3]
         
        if ch  not in chain.keys():
            t1='Error: chainID (%s) is not in coordinate. TLS group=%d !' %(ch,nt)
            perror(t1)

        if ch in chain.keys() and  n1 not in chain[ch]:
            t1='Warning: residue number (%d) is not in coordinate. TLS group=%d.' %(n1,nt)
            perror(t1)

        if ch in chain.keys() and  n2 not in chain[ch]:
            t1='Warning: residue number (%d) is not in coordinate. TLS group=%d.' %(n2,nt)
            perror(t1)

    check_residue_range(all_range)
    if id==0 : tlsxc_origin(pdb, 0)  #check TLS origin, only pdbfile and refmac

##########################################################
def check_residue_range(all_range) :
    ''' check residue range overlap :
    four type of overlaps for residue ranges (n1, n2) and (k1, k2)
    (n1, k1, k2, n2;   k1, n1, n2, k2;  n1, k1, n2, k2;   k1, n1, k2, n2;)
    '''

    nr=len(all_range)
    for i in range(0,nr):
        if len(all_range[i])<4 : continue
        n1, n2, nt1 = all_range[i][1], all_range[i][2], all_range[i][3]
       
        if n2<n1 :
            t = 'Error: Residue range problem [%d %d] in TLS group %d.' %(n1, n2, nt1)
            perror(t)
        for j in range(i+1, nr):
            if (all_range[i][0] != all_range[j][0]) : continue
            k1, k2 ,nt2 = all_range[j][1], all_range[j][2], all_range[j][3]
            ch = all_range[j][0]
            if ((n1<=k1 and n2>=k2) or     # 1
                  (k1<=n1 and k2>=n2) or     # 2
                  (k1>=n1 and n2>=k1 and k2>=n2) or # 3
                  (n1>=k1 and k2>=n1 and n2>=k2)):  # 4

                t="Error: TLS group overlaps: group=%d, residue (%d  %d) with group=%d, residue (%d  %d) ; chainID=%s " %(nt1, n1, n2, nt2, k1,k2, ch)
                perror(t)

    
##########################################################
def remove_pdb_item(fp_inp, id, item):
    '''
    id=0, fp is a pdbfile & return a pdb file;
    id=1, fp is a list of pdb & return two lists
    item, is the items to be removed.
    '''
    
    n1, n2, tls, pdb = 0, 0, [], []
    
    if id==0:
        fp=open(fp_inp, 'r').readlines()
    else:
        fp=fp_inp
        
    if item == 'tls':
        for n,x in enumerate(fp):
            if 'REMARK   3  TLS DETAILS' in x[:26] :
                n1=n+1
            elif 'S31:' in x[13:20] :
                n2=n
            elif n2>0 and 'REMARK   3' not in x[:15] :
                break
        
    pdb.extend(fp[:n1])
    tls.extend(fp[n1:n2+1])
    pdb.extend(fp[n2+1:])

    #for x in tls: print (x.strip())
    
    if id==0:
        ftmp = fp_inp + 'tmptls'
        fw = open(ftmp, 'w')
        for x in pdb: fw.write(x)
        fw.close()
        return tls, ftmp
    
    return tls, pdb
    

##########################################################
def chain_res_list(fr, id):
    '''use dic to contain chain-ID and residue range, pdb is a list.
    id=0: not include waters; id=1: include waters.
    '''
    
    chain, chain_range = {}, {}
    ch_old, n_old = '?', -99999
    
    frn=[]
    for ln in fr:
        if('ATOM' in ln[:4] or 'HETATM' in ln[:6] ): frn.append(ln)

    for ln in frn:
        if('ATOM' in ln[:4] or 'HETATM' in ln[:6] ): 
            if id==0 and 'HOH' in ln[16:20] : continue
            ch=ln[20:22].strip()
            if ch not in chain.keys() and ch != ' ' : chain[ch]=[]
            nres=ln[22:26]
            if not is_number(nres):
                perror('Warning: residue sequence (%s_%s) is not a number' %(ch,nres))
                continue
            n=int(nres)
            if ch in chain.keys():
                if(n == n_old and ch == ch_old) :
                    continue
                chain[ch].append(n)
            n_old=n
            ch_old=ch
            
    for x in chain.keys():
        if(len(chain[x])>0): chain_range[x]= [min(chain[x]), max(chain[x])]

    return chain, chain_range



##########################################################
def parse_tls_range_refmac(x, ntls):
    '''parse refmac tls range; return (ch,nres1,nres2) return none if bad!
    
    '''
    
    res=[]
    tmp = x.split(':',1)[1].strip().split()
#    print tmp,is_number(tmp[3]), int(tmp[3]), len(tmp) 
    if (len(tmp)!= 4 or tmp[0] != tmp[2] or  ' NULL ' in x  or len(str(tmp[0])) >1
        or (is_number(tmp[1]) and is_number(tmp[3]) and float(tmp[3])< float(tmp[1]))):
        
        t='Error: Wrong residue range for tls group=%s: (%s)' %(ntls, x[29:].strip())
        perror(t)
        
    elif not is_number(tmp[1]) or not is_number(tmp[3]):
        print(('Note: maybe inserted code (not treated currently): %s' % tmp))
        
    else:
        a, b=tmp[1].split('.'), tmp[3].split('.')
        res=[tmp[0], int(a[0]), int(b[0])]
        
    return res

##########################################################
def get_residue_range_from_resname(ch, res, fp):

    range=[]
    for x in fp:
        if (('ATOM' in x[:4] or 'HETA' in x[:4]) and
            (ch==x[20:22].strip() and res==x[17:20].lower())): range.append(x[22:26])
    if not range:
        perror( '\nError: the chainid/resname (%s %s) is not in the coordinate.' %(ch, res))
        return range
    
    nres1, nres2 =int(range[0]), int(range[-1])
    range1=[ch, nres1, nres2]
    
    return range1

##########################################################
def parse_tls_range_phenix_special(tls, chain, chain_range,fp):
    '''It involves in resname : for temp use. Modify later!
    '''
    
    range=[]
    for i, x in enumerate(tls):
        if pattern(i, tls, 'chain', '?', 'and', 'resname' ) :
            ch, resname=tls[2], tls[i+4]
            range_t=get_residue_range_from_resname(ch, resname, fp)
            range.append(range_t)
#            print 'Got it  %s %s ' %(ch, resname), tls
        elif pattern(i, tls, '(', 'chain', '?', 'or'):
            chs=[]
            tlstmp=tls[i:]
            for j, y in enumerate(tlstmp):
                if 'chain' in y : chs.append(tlstmp[j+1])
                if 'resname' in y and '(' not in tlstmp[j+1]:
                    resname=tlstmp[j+1]
                    for c in chs:
                        range_t=get_residue_range_from_resname(c, resname, fp)
                        range.append(range_t)
                    break
    print('Special range=%s' % range)  
    return range

##########################################################
def parse_tls_range_phenix(tlsn_in, ng, chain, chain_range, fp):
    '''tlsn: all the tls range as one string. ng: the tls group number
    chain:  chainID with a list of residue number.
    chain_range: chainID with tls residual range.
    '''
    tlsn=' '.join(tlsn_in.split())
    if len(tlsn.split('(')) != len(tlsn.split(')')) :
        t='\nError:  parentheses is not closed for TLS group = %s' %ng
        perror(t)
    
    tls, tls_range = [],[]
    s0=' '.join(tlsn.strip().replace("'", '').replace('"', '').split())
    s1=myreplace_case('resseq', 'resid', s0)
    s2=myreplace_case(' through ', ':', s1)
    tmp=s2.replace(' : ', ':').replace(': ', ':').replace(' :', ':').split()
    for i, x in enumerate(tmp):
        t=x
        if  x.endswith(':') and i<len(tmp)-1 and is_number(tmp[i+1]):
            t=x + tmp[i+1]
        if i>0 and is_number(x) and tmp[i-1].endswith(':'): continue
        
        if i>0 and 'resid' in tmp[i-1]: #check inserted code & rid of char
            t=''.join(filter(lambda y: y.isdigit() or y==':' or y=='-', x))
            
        if t and len(t) != 1 and len(t) != 2 and not is_number(t[0]) and t[0] != ':':
            tls.append(t.lower())
        else: #chid
            tls.append(t)

    
    if 'resname' in tls:
        return parse_tls_range_phenix_special(tls, chain, chain_range,fp)

    nn=len(tls)
    for i, x in enumerate(tls): # error detect
        if 'null' in x :
            t='Error: NULL values for residue range with TLS group = %s' %ng
            perror(t)
            
        elif 'chain' in x : 
            if i>nn-2 :
                t='Error: No chainID is given (tls group=%s)' %( ng)
                perror(t)
            
            else:
                chid=tls[i+1]
                if len(tls[i+1])>2 or len(tls[i+1])<1:
                    t='Error: chainID (%s) is wrong (tls group=%s)' %(tls[i+1], ng)
                    perror(t)
                else:
                    if (chid not in chain_range.keys()):
                        t='Error: chain ID (%s) not in coordinates (tls group=%s)' %(chid, ng)
                        perror(t)
                    
        elif 'segid' in x and i<nn-1:
            chid=tls[i+1]
            if len(tls[i+1])>1:
                t='Error: segidID (%s) is used (tls group=%s)' %(tls[i+1], ng)
                perror(t)
            
    for i, x in enumerate(tls):
        if 'all' in x:
            for y in chain_range.keys():
                tls_range.append([y,chain_range[y][0], chain_range[y][1]])
            
        elif pattern(i, tls, 'chain', '?', 'and', 'not' , 'resid') :
            ch=tls[i+1]
            if ch  in chain.keys():
                tt = parse_tls_phenix_not(i,ch, tls, chain[ch], chain_range, 0)
                tls_range.extend(tt)
                
        elif pattern(i, tls, 'chain', '?', 'and', 'not' , '(', 'resid') :
            ch=tls[i+1]
            if ch  in chain.keys():
                tt = parse_tls_phenix_not(i, ch,tls, chain[ch], chain_range, 1)
                tls_range.extend(tt)
            
        elif 'chain' in x and i<nn-1:
            chid=tls[i+1]
            
            if (chid not in chain_range.keys()):
                continue
            
#            elif (pattern(i, tls, 'chain', '?') and (nn==2 or nn==4)) : #one chain
            elif (pattern(i, tls, 'chain', '?', 'or') or
                  pattern(i, tls, 'chain', '?',')', 'or') or
                  (pattern(i, tls, 'chain', '?') and i>nn-3 ) or
                  pattern(i, tls, 'chain', '?', ')') ) : #one chain
                
                n1, n2 = chain_range[chid][0], chain_range[chid][1]
                tls_range.append([chid, n1,n2])
                
                
            elif (pattern(i, tls, 'chain', '?', 'or', 'chain', '?') ): #two chain
                
                n1, n2 = chain_range[chid][0], chain_range[chid][1]
                tls_range.append([chid, n1,n2])

                ch2=tls[i+4]
                if (ch2 not in chain_range.keys()): continue
                n1, n2 = chain_range[ch2][0], chain_range[ch2][1]
                tls_range.append([ch2, n1,n2])

              
            elif (i<nn-4 and (pattern(i, tls, 'chain', '?', 'and', 'res'))):
                
                n1, n2 = get_range_ph(chain_range[chid], tls[i+4])
                tls_range.append([chid, n1,n2])
                
            elif ((pattern(i, tls, 'chain', '?', 'and', '(', 'res'))):
                get_range_phe(ng, chid, i+4, 1, chain_range, tls, tls_range)

            elif ((pattern(i, tls, 'chain', '?', 'and', '(', '(', 'res'))):
                get_range_phe(ng, chid, i+5, 2, chain_range, tls, tls_range)


        elif i==0 and  pattern(i, tls, 'res', '?', 'and', 'chain', '?') :
            chid=tls[i+4]
            if (chid not in chain_range.keys()): continue
            n1, n2 = get_range_ph(chain_range[chid], tls[i+1])
            tls_range.append([chid, n1,n2])
            
        elif i==0 and  pattern(i, tls, '(', 'res', '?',')', 'and', 'chain', '?') :
            chid=tls[i+6]
            if (chid not in chain_range.keys()): continue
            
            n1, n2 = get_range_ph(chain_range[chid], tls[i+2])
            tls_range.append([chid, n1,n2])

        elif (i==0 and pattern(i, tls, '(','res', '?', 'or', 'res', '?')
              and pattern(8, tls, 'chain',  '?' )) :
            chid=tls[tls.index('chain') +1]
            if (chid not in chain_range.keys()): continue
            get_range_phe(ng, chid, i+0, 1, chain_range, tls, tls_range)

        elif (pattern(i, tls, 'res') and 'chain' not in tls ):
            if len(chain_range.keys())>1:
                t='Error: wrong residue selection, no chain identifer (tls group=%s).' %(ng)
                perror(t)
            else:
                chid=list(chain_range.keys())[0]
                n1, n2 = get_range_ph(chain_range[chid], tls[i+1])
                tls_range.append([chid, n1,n2])
                
    check_parsed_tls_ranges(tls_range, chain, ng)                    
        
    return tls_range                                    
       
##########################################################
def check_parsed_tls_ranges(tls_range, chain, ng):
    
    if not len(tls_range):
        t='Error: failed to parse residue range (TLS group=%s).' %ng
        perror(t)
    else: #check if residue exist
        for x in tls_range:
            if not x : continue
            k1, v1, v2 = x[0], x[1], x[2]
            if k1 not in chain.keys():
                t='Error: chain ID (%s) not in coordinates (tls group=%s)' %(k1, ng)
                perror(t)
            else:
                if v1 not in chain[k1]:
                    t='Warning: residue (%s_%d) not in coordinates (tls group=%s)' %(k1,v1, ng)
                    perror(t)
                if v1 != v2 and v2 not in chain[k1]:
                    t='Warning: residue  (%s_%d) not in coordinates (tls group=%s)' %(k1,v2, ng)
                    perror(t)
    
##########################################################
def tls_selection_one_string(fp, i, nc):
    ''' fp: a list;  i: the position for SELECTION; nc, start column of range.
    return a string for this selection.
    '''
    
    tls = fp[i].split(':',1)[1].strip() + ' '  # the first one
#    tlsn=tls.strip().replace('(', ' ( ').replace(')', ' ) ')
    tlsn=''
    ss=''
    for n  in range(i+1, len(fp)):
        if 'REMARK   3    ORIGIN FOR ' in  fp[n][:25]:  #stop here
            tls = tls + ss
            tlsn=tls.strip().replace('(', ' ( ').replace(')', ' ) ')
            break
        ss = ss + (fp[n][nc:].strip()  + ' ')

    return tlsn

##########################################################
def pattern(i, tls, *ss):
    '''the pattern of TLS after position i.  ss is a list of input
    '''
    
    key, nn = 1, len(tls)
    
    if len(ss)>nn-i : return 0
    for k, x in enumerate(ss):
        n=i+k
        if n > nn-1 or x == '?' or (x=='res' and x==tls[n][:3]) : continue
        if x != tls[n]:
            key=0
            break
        
    return key
    
##########################################################
def parse_tls_phenix_not(i,ch, tls, chain, chain_range, id):
    '''parse tls range involve NOT;
    id=0 for style 'chain F and not resid 0'
    id=1 for style 'chain F and not  not (resseq 539:581 or ..)'
    '''
    
    chain_range_not=[]
    chain_range_yes=[]
    chain_new=chain
    
    if id==0:
        n1, n2= get_range_ph(chain_range, tls[i+5])
        chain_range_not.append([n1, n2])
        
    else:
        m, tmp_list = get_parenthetic_contents(i+4, tls)
        for i, x in enumerate(tmp_list):
            if 'resid' in x :
                n1, n2= get_range_ph(chain_range, tmp_list[i+1])
                chain_range_not.append([n1, n2])

    chain_range_not.sort()
    
    nn1=len(chain_range_not)
    chain_new.sort()
    nn2=len(chain_new)
    for i, x in  enumerate(chain_range_not):
        if i ==0:  #first one 
            if x[0] - chain_new[0] >2:
                chain_range_yes.append([ch, chain_new[0],  x[0]-1])
        if i==nn1-1: #last one
            if chain_new[nn2-1] - x[1] >2:
                chain_range_yes.append([ch, x[1]+1, chain_new[nn2-1]] )

        if i<nn1-1 and chain_range_not[i+1][0] - x[1]>2 :
            chain_range_yes.append([ch, x[1]+1, chain_range_not[i+1][0]-1 ])
            
#    print('chain_range_not=',chain_range_not)
#    print('chain_range_yes=',chain_range_yes)
    return chain_range_yes
            
    
##########################################################
def get_parenthetic_contents(start, lists):
    ''' start:  the start position;  lists: input a list or string
    return the position of last ')' and the contents in ()
    '''
    
    nn=len(lists)
    tmp_list=[]
    m1, m2, n=0, 0, 0
    for i in range(start, nn):
        if lists[i] == '(' : m1=m1+1
        if lists[i] == ')' : m2=m2+1
        if type(lists) == type(list()): # a list
            tmp_list.append(lists[i])
            
        if (m1>0 and m1 == m2) :
            if type(lists) == type(list()): # a list
                return n, tmp_list
            elif type(lists) == type(str()): # a str
                return n, lists[start:n]
            break
        n=n+1
       
    if m1 != m2:
        ss=' '.join(lists)
        t= 'Error: Parenthese not equal in line "%s" \n' %ss
        perror(t)
    return 0, tmp_list
    

##########################################################
def get_range_phe(ng, chid, st, nb, chain_range, tls, tls_range):
    '''st: start point;  nb: number of bracket;  chain_range: {ch:[n1,n2]}
    tls_range: return the value 
    '''
    
    nn=len(tls)
    kk=0
    for m in range(st, nn):
        if (nb==1 and tls[m] == ')') or (nb==2 and tls[m] == ')'and m+1<nn and tls[m+1] == ')' ):
            kk=1
            break
        if ':' in tls[m]:
            n1, n2 = get_range_ph(chain_range[chid], tls[m])
            tls_range.append([chid, n1,n2])
            
    if kk==0:
        ss = 'Error: residue range is possibly truncated (TLS group=%s)!' %ng
        perror(ss)
    
##########################################################
def get_range_ph(chain_range, tls):
        
    n1, n2 = -9999, -9999
    if ':' in tls:
        if len(tls.split(':'))==2:
            t=tls.split(':')
            if is_number(t[0]) : n1 = int(t[0])
            if is_number(t[1]) : n2 = int(t[1])
    else:
        if '-' in tls and len(tls.split('-'))==2:
            t=tls.split('-')
            if is_number(t[0]) : n1 = int(t[0])
            if is_number(t[1]) : n2 = int(t[1])
        else:
            if is_number(tls): n1=n2=int(tls)

    if n1== -9999 : n1=chain_range[0]
    if n2== -9999 : n2=chain_range[1]

    return n1, n2

##########################################################
def parse_tls_range_buster(tls, ng,  chain, chain_range):
    '''tls: all the tls range as one string.  ng: the tls group number
    chain:  chainID with a list of residue number.
    chain_range: chainID with tls residual range.
    
    '''

    tls_range = []

    tlsn=tls.replace(' |', '|').replace('| ', '|').replace('{','').replace('}','')   
    
    if '|' not in tlsn or 'null' in tlsn.lower():
        t='Error: problem with TLS group=%s. \n' %(ng)
        perror(t)

    elif ' - ' not in tlsn : # one residue or all in one chain
        st=tlsn.split()
        
        for x in st:
            if '|' in x :
                tmp=parse_tls_range_buster_1(x,  chain_range)
                if x[0]=='*' :
                    tls_range = tmp
                else:
                    tls_range.append(tmp)
                                    
    else : # multiple ranges
        ss=tlsn.split()
        i=0
        nn=len(ss)
        for k in range(i,nn):
            if '|' in ss[k] and '*'  in ss[k] and nn-k>1 and '-' not in ss[k+1] : #isolate entity
                tmp=parse_tls_range_buster_1(ss[k], chain_range)
                if tmp :
                    tls_range.append(tmp)
                else:
                    t='Error: residue range (%s) not parsed (tls group=%s)' %(ss[k], ng)
                    perror(t)
                    
            elif '|' in ss[k] and nn-k>1 and '-' not in ss[k+1] : #isolate entity
                tmp=parse_tls_range_buster_1(ss[k], chain_range)
                
            elif  '|' in ss[k] and nn-k>2 and '-' in ss[k+1] :
                
                ch1, res1 = parse_tls_range_buster_2(ss[k])
                
                if '|' not in ss[k+2] :
                    if is_number(res1) and is_number(ss[k+2]):
                        tls_range.append([ch1, int(res1), int(ss[k+2])])
                    else:
                        t='Error: residue range (%s - %s) not parsed (tls group=%s)' %(ss[k],ss[k+2], ng)
                        perror(t)
                        
                else:
                    ch2, res2 = parse_tls_range_buster_2(ss[k+2])
                    if (ch1 == ch2 and is_number(res1) and is_number(res2) ) :
                        tls_range.append([ch1, int(res1), int(res2)])
                    else:
                        print('Warning: parse problem for TLS group=%s' %ng)

                i=i+2

            i=i+1
                
    check_parsed_tls_ranges(tls_range, chain, ng)
    
    return tls_range

##########################################################
def parse_tls_range_buster_2(x):

    s = x.split('|')
    ch, res = s[0], s[1]
    return ch, res
                    
##########################################################
def parse_tls_range_buster_1(x, chain_range):
    ''' parse single items  (A|*,  A|45, A|12-45, A|12-A|45)
    '''

    tmp=[]
    ch, res = parse_tls_range_buster_2(x)
    if '*' in ch  and '*' in res: # all chains
        for k, v in chain_range.items():
            tmp.append([k, v[0], v[1]])
        return tmp
    
    if ch not in chain_range.keys() :
        t='Error: chain ID (%s) not in coordinates ' %(ch)
        perror(t)
        return tmp
    
    if '*' in x:  
        tmp=[ch, chain_range[ch][0], chain_range[ch][1]]
        
    elif  is_number(res):
        tmp=[ch, int(res),int(res)]
        
    else:  #for connected case (A|12-45, A|12-A|45)
        tt=x.split('|')
        if '-' in res and len(tt)==2 : #A|12-45
            t=res.split('-')
            if len(t) !=2 : return tmp
            t1, t2 = t[0], t[1]
            if is_number(t1) and is_number(t2): tmp=[ch, int(t1),int(t2)]
            
        elif '-' in res and len(tt)==3:  #A|12-A|45
            t1=x.split('-')
            if len(t1)!=2 : return tmp
            
            if '|' in t1[0] and '|' in t1[1]:
                res1, res2 = t1[0].split('|'), t1[1].split('|')
                
                if is_number(res1[1]) and is_number(res2[1]) and res1[0] == res2[0]:
                    tmp=[ch, int(res1[1]),int(res2[1])]
                else:
                    t='Error: The first chain ID (%s) differs from the second(%s)' %(res1[0],x)
                    perror(t)
        
    return tmp

##########################################################
def tlsxc_origin(pdbfile, info):
    ''' calculate the origin for TLS groups (occupancy weighted)
    refmac & buster: center of coordinate.
    phenix: center of mass
    info ==1, print all. info==0, only print error/warning
    
    '''

    fp=open(pdbfile, 'r').readlines()

    pdbid, prog, tls, xyz = '', '?', 0, [] 
    for x in fp:
        if 'HEADER' in x[:6]:
            pdbid=x[62:66].strip()
            
        elif 'REMARK   3   PROGRAM     :' in x:
            prog=x.split(':')[1].strip()
            
        elif 'S21:' in x and 'S22:' in x:
            tls = tls + 1
            
        elif 'ATOM  ' in x[:6] or 'HETATM' in x[:6]:
            xyz.append(x)
            
        elif 'ENDMDL' in x[:6] :
            print('Warning. This is a multiple model entry (first model used!).\n')
            break

    if tls==0:
        print('Note: No TLS groups are detected in file=(%s).' %pdbfile)
        return
    
    chain, ch_range = chain_res_list(fp, 1)
    nxyz, ntls, nph, nref, nbu = [], '?', 0, 0, 0
    for i, x in enumerate(fp):
        if 'REMARK' in x[:6] and 'BUSTER' in x[24:35]:
            nbu=1
            
        elif 'TLS GROUP :' in x [:25]:
            ntls=x.split(':')[1].strip()
            
        elif 'REMARK   3    RESIDUE RANGE :' in x [:29]:  #refmac
            nref = nref + 1
            if nref == 1 and info==1:
                print('\nNOTE: TLS format is REFMAC. Origin is OCC weighted center of coordinate.')
                
            t = parse_tls_range_refmac(x, ntls)
            if not t: continue
            tmp = [t[0], t[1], t[0], t[2]]
            
            txyz = get_xyz_tls(tmp, xyz)
            nxyz.extend(txyz)
            
        elif 'REMARK   3    SELECTION:' in x [:29] and 'phenix' in prog.lower():  #phenix,
            nph = nph + 1
            if nph == 1 and info==1:
                print('\nNOTE: TLS format is PHENIX. Origin is OCC weighted center of mass.')
            
            tls_str = tls_selection_one_string(fp, i, 25)
            if not tls_str : continue
            t = parse_tls_range_phenix(tls_str, ntls, chain, ch_range,fp)
            
            for z in t:
                if len(z) <3 : continue
                tmp=[z[0], z[1], z[0], z[2]]
                txyz=get_xyz_tls(tmp, xyz)
                nxyz.extend(txyz)
            
        elif (('REMARK   3    SET :' in x [:20]) or
              ('REMARK   3    SELECTION:' in x[:29] and 'buster' in prog.lower() ))  :  #buster
            nbu = nbu + 1
            if nbu == 1 and info==1:
                print('\nNOTE: TLS format is BUSTER. Origin is OCC weighted center of coordinate.')
                print('        The origin column size is not the same as phenix & refmac. ')
                
            if 'SELECTION:' in x  and '{' not in x[20:]: #buster (tls in phenix/refmac format)
                t = parse_tls_range_refmac(x, ntls)
                if not t: continue
                tmp = [t[0], t[1], t[0], t[2]]
                txyz = get_xyz_tls(tmp, xyz)
                nxyz.extend(txyz)
                
            else:
                if '{' in  x[20:]:
                    tls_str = tls_selection_one_string(fp, i, 25)
                else:
                    tls_str = tls_selection_one_string(fp, i, 20)
                    
                t = parse_tls_range_buster(tls_str, ntls, chain, ch_range)
                for z in t:
                    if len(z) <3 : continue
                    tmp=[z[0], z[1], z[0], z[2]]
                    txyz=get_xyz_tls(tmp, xyz)
                    nxyz.extend(txyz)
            
        elif 'REMARK   3    ORIGIN FOR THE GROUP' in x:
            xc=x.rstrip().split(':')[1].split()
      
            if len(xc)!=3:                 
                xc=[x[39:48], x[48:57], x[57:66]]
                if nbu >0 : xc=[x[39:49], x[49:59], x[59:69]]

            if (not is_number(xc[0]) or not is_number(xc[1]) or not is_number(xc[2])):
                t='Error: Wrong TLS origin (%s) for group=%s' %(xc,ntls)
                perror(t)
                continue
            
            nxc=[]
            oxc=[float(y) for y in xc]
            if nref:
                nxc=calc_tlsxc(nxyz, 'refmac')
            elif nph:
                nxc=calc_tlsxc(nxyz, 'phenix')
            elif nbu:
                nxc=calc_tlsxc(nxyz, 'buster')

          #  print xc, nxc
            if(len(nxc)<1):
                t='\nWarning: Origin not calculated. No xyz found in TLS group=%s' %ntls
                perror(t)
            else:
                dif=[math.fabs(nxc[0]-oxc[0]), math.fabs(nxc[1]-oxc[1]), math.fabs(nxc[2]-oxc[2])]
                if info==1:
                    print('\nTLS origin (calculated): %8.3f %8.3f %8.3f' %(nxc[0], nxc[1],nxc[2]))
                    print('TLS origin (reported  ): %8.3f %8.3f %8.3f' %(oxc[0], oxc[1],oxc[2]))
                if( max(dif)>2 ):
                    t="Warning: (%s:) TLS origin differs from calculated(max_dev=%.2f),group=%s." %(pdbid,max(dif), ntls)
                    perror(t)
                else:
                    if info==1:
                        print("(%s:) TLS origin is similar to the calculated for group=%s." %(pdbid,ntls))
                
            nxyz=[]
            
        elif 'MODEL' in x[:6]:
            print('Warning. This is a multiple model entry!\n')
        elif 'ATOM' in x[:6]:
            break
    #for x in xyz : print(x)

##########################################################
def calc_tlsxc(nxyz, id):
    ''' get the center of coordinate (refmac) and center of mass (phenix)
    
    '''

    occup=0.0
    for x in nxyz:
        if len(x.strip()) > 60: occup = occup + float(x[55:60])
        
    if occup <=0:
        return []
    
    occmass, xt, yt, zt = 0.0, 0.0, 0.0, 0.0
    for ln in nxyz:
        x, y, z = float(ln[29:38]), float(ln[38:46]), float(ln[46:54])
        occ = float(ln[55:60])
        mass=1
        if id.lower()=='phenix' :
            atom_mass = 0.0
            if len(ln)>78 :  atom_mass = mass_of_atom(ln[76:78].strip())
            if atom_mass > 0 : mass = atom_mass
            
        xt=xt + x*occ*mass
        yt=yt + y*occ*mass
        zt=zt + z*occ*mass
        occmass = occmass + occ*mass

    return [xt/occmass, yt/occmass, zt/occmass]
    
##########################################################
        
def get_xyz_tls(res, xyz):
    '''input res (ch, range1, ch, range2); xyz coordinates.
    '''
    nxyz=[]
    for x in xyz[:]:
        if not atom_record(x) : continue
        if (res[0] == x[20:22].strip() and is_number(x[22:26]) and  int(res[1]) <= int(x[22:26]) <= int(res[3])):
            nxyz.append(x)
#            xyz.remove(x)
    return nxyz


##########################################################
def atom_record(x):
    val=0
    if ('ATOM' in x[:4] or 'HETA' in x[:4]  or 'ANISOU' in x[:6]): val=1
    
    return val

##########################################################
def residue(resname):

    val=0
    resid=['GLY','ALA','VAL','LEU','ILE','CYS','MET','PRO','PHE','TRP','TYR',
           'HIS','LYS','ARG','ASP','GLU','ASN','GLN','THR','SER','MSE','UNK',
           '  A','  G','  C','  T',' DA',' DG', ' DC',' DT','  U','  I']

    if resname in resid: val=1
    return val

##########################################################
def mass_of_atom(atom):
    atom_mass={
        "H": 1.0079 ,
        "HE": 4.0026 ,
        "LI": 6.941 ,
        "BE": 9.01218 ,
        "B": 10.81 ,
        "C": 12.011 ,
        "N": 14.0067 ,
        "O": 15.9994 ,
        "F": 18.9984 ,
        "NE": 20.179 ,
        "NA": 22.98977 ,
        "MG": 24.305 ,
        "AL": 26.98154 ,
        "SI": 28.0855 ,
        "P": 30.97376 ,
        "S": 32.06 ,
        "CL": 35.453 ,
        "AR": 39.948 ,
        "K": 39.0983 ,
        "CA": 40.08 ,
        "SC": 44.9559 ,
        "TI": 47.9 ,
        "V": 50.9415 ,
        "CR": 51.996 ,
        "MN": 54.938 ,
        "FE": 55.847 ,
        "CO": 58.9332 ,
        "NI": 58.7 ,
        "CU": 63.546 ,
        "ZN": 65.38 ,
        "GA": 69.72 ,
        "GE": 72.59 ,
        "AS": 74.9216 ,
        "SE": 78.96 ,
        "BR": 79.904 ,
        "KR": 83.8 ,
        "RB": 85.4678 ,
        "SR": 87.62 ,
        "Y": 88.9059 ,
        "ZR": 91.22 ,
        "NB": 92.9064 ,
        "MO": 95.94 ,
        "TC": 98.0 ,
        "RU": 101.07 ,
        "RH": 102.9055 ,
        "PD": 106.4 ,
        "AG": 107.868 ,
        "CD": 112.41 ,
        "IN": 114.82 ,
        "SN": 118.69 ,
        "SB": 121.75 ,
        "TE": 127.6 ,
        "I": 126.9045 ,
        "XE": 131.3 ,
        "CS": 132.9054 ,
        "BA": 137.33 ,
        "LA": 138.9055 ,
        "CE": 140.12 ,
        "PR": 140.9077 ,
        "ND": 144.24 ,
        "PM": 145.0 ,
        "SM": 150.4 ,
        "EU": 151.96 ,
        "GD": 157.25 ,
        "TB": 158.9254 ,
        "DY": 162.5 ,
        "HO": 164.9304 ,
        "ER": 167.26 ,
        "TM": 168.9342 ,
        "YB": 173.04 ,
        "LU": 174.967 ,
        "HF": 178.49 ,
        "TA": 180.9479 ,
        "W": 183.85 ,
        "RE": 186.207 ,
        "OS": 190.2 ,
        "IR": 192.22 ,
        "PT": 195.09 ,
        "AU": 196.9665 ,
        "HG": 200.59 ,
        "TL": 204.37 ,
        "PB": 207.2 ,
        "BI": 208.9804 ,
        "PO": 209.0 ,
        "AT": 210.0 ,
        "RN": 222.0 ,
        "FR": 223.0197 ,
        "RA": 226.0254 ,
        "AC": 227.0278 ,
        "TH": 232.0381 ,
        "PA": 231.0359 ,
        "U": 238.029 ,
        "NP": 237.0482 ,
        "PU": 244.0 ,
        "AM": 243.0 ,
        "CM": 247.0 ,
        "BK": 247.0 ,
        "CF": 251.0 ,
        "ES": 252.0 ,
        "FM": 257.0 ,
        "MD": 258.0 ,
        "NO": 259.0 ,
        "LR": 260.0 ,
        "RF": 261.0 ,
        "DB": 262.0 ,
        "SG": 263.0 ,
        "BH": 262.0 ,
        "HS": 265.0 ,
        "MT": 266.0 ,
        "CN": 277.0 , 
        "X": 12.011 
        
        }
    #for x in atom_mass.keys(): print(x, atom_mass[x])
    
    if atom.upper() in atom_mass.keys():
        return atom_mass[atom.upper()]
    else:
        print('Warning: The atom (%s) not exist in the list' %atom)
        return 0
    


##########################################################
# Below is everything about correcting B factors using TLSANL
##########################################################
##########################################################
def proc_tlsanl(pdbfile, outfile):
    '''run tlsanl program: correct residue ranges if they are wrong.
    '''
    
    delete_file(outfile)
    print ('%s' %(60*'-'))
    (pdb_tls1, log1, pdb_tls2, log2, tls_inp, tls_new_inp, pdb_0tls,
     pdb_new, pdb_newb) = ('TMP.LOG',)*9
                
    tls = precheck_tls(pdbfile) 
    if(tls == 2): # phenix
        print('Note:  TLS format is phenix (%s).' %pdbfile)
        return 0
    elif (tls == 3): # Buster
        print('Note:  TLS format is Buster (%s).' %pdbfile) 
        return 0

    tls_inp = get_tlsinp(pdbfile)
    pdb_tls1, log1 ='',''
    if tls==1: (pdb_tls1, log1) = do_tlsanl(pdbfile, tls_inp, 't', 'full')
    
    if os.path.exists(pdb_tls1) and os.path.getsize(pdb_tls1)>500:
        shutil.move(pdb_tls1, outfile)
        print ("\n%s; TLSANL was successful (use original residue-range).\n" %pdbfile)
    else:
        pdb_0tls = delete_0_origin(pdbfile) #remove tls group with 0 origin
        pdb_new = correct_tls(pdb_0tls)  #do various corrections
        tls_new_inp = get_tlsinp(pdb_new) #extract TLS groups into tlsanl
        (pdb_tls2, log2) = do_tlsanl(pdb_new, tls_new_inp, 't', 'full')
        
        if os.path.exists(pdb_tls2) and os.path.getsize(pdb_tls2)>500:
            shutil.move(pdb_tls2, outfile)
            print ("\n%s; TLSANL was successful (residue-range corrected).\n" %pdbfile)
        else:
            t="Error: TLS correction failed (%s).\n" %pdbfile
            perror(t)
            if  os.path.exists(log2) and os.path.getsize(log2)>100 : 
                fp=open(log2, "r").readlines()
                n=0
                for ln in fp:
                    n=n+1
                    if "ERROR:" in ln:
                        print(pdbfile + '; ' + ln)
                        if len(fp[n].strip())>2:  print(pdbfile + '; ' + fp[n])
                        break
            
    wk=0    
    if os.path.exists(outfile) and os.path.getsize(outfile)>100 :
        #bigb, bval = check_bfactor(outfile) #only show Biso status
        #finalize_pdb(pdbfile, outfile)
        #print('The pdbfile after tlsanl correction = %s\n' %outfile)
        wk=1
        #if bigb>0: 
        #    wk=0
            #delete_file(outfile)
            
    delete_file(pdb_tls1, log1, pdb_tls2, log2, pdb_0tls, pdb_new, pdb_newb)
    delete_file(tls_inp, tls_new_inp)
    print ('%s' %(60*'-'))
                
    return wk  

##########################################################
def precheck_tls(pdbfile):
    '''pre-check the type of TLS and possible errors
    '''
    
    wk=1   #default refmac
    
    fr=open(pdbfile, "r")
    for x in fr:
            
        if "REMARK   3    RESIDUE RANGE : " in x:
            tmp=x.split(':')[1].split()
            if(len(tmp) !=4 or is_digit(tmp[1])==0 or  is_digit(tmp[3])==0 or
               (int(tmp[1])<0 and (int(tmp[3])==-1 or int(tmp[3])>=999))):
                wk=0
                break
            
        elif "REMARK" in x[:6] and  " SELECTION: " in x[13:25] :
            wk=2
            break
            
        elif "REMARK" in x[:6] and  " SET : " in x[12:20] :
            wk=3
            break
            
        elif 'ATOM' in x[:4] or 'HETATM' in x[:6]:
            break
        
    fr.close()

    
    return wk

##########################################################
def get_tlsinp(pdbfile):
    '''Get TLS input for TLSANL:
    This could be done by tlsextract, but data columns can be connected,
    then TLSANL will fail.
    '''
    
    tlsinp=pdbfile + '_tls.inp'
    fr=open(pdbfile, 'r')
    fo=open(tlsinp, 'w')
    fo.write("REFMAC \n\n")

    n=0
    for ln in fr:
        if 'BULK SOLVENT MODELLING' in ln or ln[0:4]=='ATOM' or ln[0:6]=='HETATM' : break
        if 'REMARK   3' not in ln[:12] : continue
        
        if 'TLS GROUP :' in ln:
            ntls = get_value_after_id(ln,":")
            fo.write('TLS GROUP : ' + ntls + '\n')
            
        elif 'RESIDUE RANGE :' in ln: #for refmac
            t1=parse_tls_range_refmac(ln, ntls)
            if len(t1)<2 : continue
            fo.write("RANGE  '%s%4d.' '%s%4d.' ALL\n" %(t1[0],t1[1],t1[0],t1[2]) )
            
        elif '  SELECTION: ' in ln: #for phenix
            if 'NULL' in ln:
                print(pdbfile + '; Error: wrong residue range(%s)\n' %ln[25:].strip())
                break
            n=n+1
            if n==1: chain, chain_range, chain1, chain_range1, chain2, chain_range2 =chain_res_range(pdbfile)
            range1=convert_range_phenix_to_refmac(ln,chain_range)
            for k, v in range1.items():
                for x in v:
                    arg="RANGE  '%s%4s.' '%s%4s.' ALL\n" %(k,x[0],k,x[1])
                    fo.write(arg)
                
            
        elif 'REMARK   3    ORIGIN FOR THE GROUP (' in ln and '):' in ln:
            tmp=ln.split()
            
            if 'NULL' in ln or  '0.0000   0.0000   0.0000' in ln or len(tmp)>10:
                print(pdbfile + '; Error: wrong TLS origin (%s).\n' %ln[40:].strip() )
                break
            
            fo.write('ORIGIN  %8s %s %8s \n' %(ln[39:48],ln[48:57],ln[57:66]) )
            
        elif 'T11:' in ln and 'T22:' in ln:
            #T11, T22 = ln[20:29], ln[34:43]
            T11 = get_one_string_after_id(ln, 'T11:')
            T22 = get_one_string_after_id(ln, 'T22:')
            
        elif 'T33:' in ln and 'T12:' in ln:
            #T33, T12 = ln[20:29], ln[34:43]
            T33 = get_one_string_after_id(ln, 'T33:')
            T12 = get_one_string_after_id(ln, 'T12:')
            
        elif 'T13:' in ln and 'T23:' in ln:
            #T13, T23 = ln[20:29], ln[34:43]
            T13 = get_one_string_after_id(ln, 'T13:')
            T23 = get_one_string_after_id(ln, 'T23:')

            if 'NULL' in T11+ T22+ T33+ T12+ T13+ T23:
                t= 'Error: "NULL" values in matrix Tij for tls group=%s\n' %ntls
                perror(t)
                break
              
            fo.write('T  %s %s %s %s %s %s\n' %(T11, T22, T33, T12, T13, T23))
            
        
        elif 'L11:' in ln and 'L22:' in ln:
            #L11, L22 = ln[20:29], ln[34:43]
            L11 = get_one_string_after_id(ln, 'L11:')
            L22 = get_one_string_after_id(ln, 'L22:')
            
        elif 'L33:' in ln and 'L12:' in ln:
            #L33, L12 = ln[20:29], ln[34:43]
            L33 = get_one_string_after_id(ln, 'L33:')
            L12 = get_one_string_after_id(ln, 'L12:')
            
        elif 'L13:' in ln and 'L23:' in ln:
            #L13, L23 = ln[20:29], ln[34:43]
            L13 = get_one_string_after_id(ln, 'L13:')
            L23 = get_one_string_after_id(ln, 'L23:')

            if 'NULL' in L11+ L22+ L33+ L12+ L13+ L23:
                t= 'Error: "NULL" values in matrix Lij for tls group=%s\n' %ntls
                perror(t)
                
                break
            
            fo.write('L  %s %s %s %s %s %s\n' %(L11, L22, L33, L12, L13, L23))

        
        elif 'S11:' in ln and 'S12:' in ln  and 'S13:' in ln :
            #S11, S12, S13 = ln[20:29], ln[34:43], ln[48:57]
            S11 = get_one_string_after_id(ln, 'S11:')
            S12 = get_one_string_after_id(ln, 'S12:')
            S13 = get_one_string_after_id(ln, 'S13:')
            
        elif 'S21:' in ln and 'S22:' in ln  and 'S23:' in ln :
            #S21, S22, S23 = ln[20:29], ln[34:43], ln[48:57]
            S21 = get_one_string_after_id(ln, 'S21:')
            S22 = get_one_string_after_id(ln, 'S22:')
            S23 = get_one_string_after_id(ln, 'S23:')
            
        elif 'S31:' in ln and 'S32:' in ln  and 'S33:' in ln :
            #S31, S32, S33 = ln[20:29], ln[34:43], ln[48:57]
            S31 = get_one_string_after_id(ln, 'S31:')
            S32 = get_one_string_after_id(ln, 'S32:')
            S33 = get_one_string_after_id(ln, 'S33:')

            if 'NULL' in S11+ S22+ S33+ S12+ S13+ S23:
                t= 'Error: "NULL" values in matrix Sij for tls group=%s\n' %ntls
                perror(t)
                break
            
            S22_S11 = float(S22) - float(S11)
            S11_S33 = float(S11) - float(S33)
            
            fo.write('S  %8.4f %8.4f %s %s %s %s %s %s\n\n'
                     %(S22_S11, S11_S33, S12, S13, S23, S21, S31, S32 ))
            

    fr.close()
    fo.close()
    return tlsinp
##########################################################
def get_one_string_after_id(line, id):
    """ get one value after a given id """
    
    if id not in line : return "?"
    value=line.split(id)[1].split()[0].strip()
    if not value : value= "?"
    #print(line, id, value)
    return value


##########################################################
def convert_range_phenix_to_refmac(ln,chain_range):
    ''' Convert TLS residue range from phenix to refmac

    '''

    line=ln[24:].replace('(', '').replace(')', '').replace("'", '') 
    tmp = line.upper().split()
    n=len(tmp)
    
    val={}
    v=[]

    if n ==2:
        if tmp[0] == 'CHAIN' and len(tmp[1])==1:
            
            v.append(chain_range[tmp[1]])
            val[tmp[1]] = v
            
    elif n==5:
        if ('CHAIN' in tmp[0] and len(tmp[1])==1 and
            'AND' in tmp[2] and 'RES' in tmp[3]):

            res = resi_range(tmp[4])
            v.append(res)
            val[tmp[1]] = v
        
        elif ('RES' in tmp[0] and len(tmp[4])==1 and
            'AND' in tmp[2] and 'CHAIN' in tmp[3]):

            val[tmp[4]] = v
            if ':' in tmp[1]:
                t1=tmp[1].split(':')
                if not t1[0]: t1[0] = '-10'
                v.append([t1[0], t1[1]])
    
    elif n==8:
        if('CHAIN' in tmp[0] and 'AND' in tmp[2] and 'RES' in tmp[3] and
           'OR' in tmp[5] and 'RES' in tmp[6]):
            chain=tmp[1]
            val[chain] = v
            res = resi_range(tmp[4])
            v.append(res)
            res = resi_range(tmp[7])
            v.append(res)
            
          
#    print(tmp , n, val)
    
    return val


##########################################################
def resi_range(tmp):
    
    v=[]
    if ':' in tmp:
        t1=tmp.split(':')
        if not t1[0]: t1[0] = '-10'
        v= [t1[0], t1[1]]
                
    elif ':' not in tmp and '-' in tmp:
        t1=tmp.split('-')
        v= [t1[0], t1[1]]
                
    else:
        v=[tmp, tmp]

    return v
    
##########################################################

def do_tlsanl(pdbfile, tls_new_inp, bresid, isoout):
    """ execute csh script to get a new pdb file with modified B=Bres+Btls """
    

    csh_script="""#!/bin/tcsh  -f

############################################################### 
#   tlsanl: to correct B factors (B_total = B_residual + B_tls),
#     when the following conditions are satisfied
#     a). The PDB coordates were refined by REFMAC
#     b). When TLS was involved in refinement
#     c). When the residue ranges are correct (so that tlsanl works).
##################### Usage ################################## 
#  tlsanl_script.csh  pdbfile      (use tls.inp  by tlsextract)
#    or 
#  tlsanl_script.csh  pdbfile  tls.inp  (use the provided tls.inp)
##################### Do TLSAN ###############################

if( $#argv <1 ) then
    echo "Usage: tlsanl_script.csh  pdbfile"
    echo " OR "
    echo "Usage: tlsanl_script.csh  pdbfile  tls.inp"
    exit()
endif

set pdbfile=$1
set refmac=`grep "^REMARK   3   PROGRAM " -m 1 $pdbfile | awk '{print $5}' `  
if( $refmac !~ "REFMAC" ) then
echo "Warning: Entry ($1) not refined by REFMAC!"
#    exit()
endif

set aniso=`grep "^ANISOU" -m 2 $pdbfile | wc -l  | awk '{print $1}' ` 
if ( $aniso >1 ) then
echo "Warning: ANISOU records exist. B_total will be overcorrected!"
#echo "Warning: ANISOU records exist. no tls correction was applied"
#    exit()
endif

if ($#argv <2) then    # no tls.inp, extract tls from PDB
    
    set tls_inp = "${pdbfile}_tls.inp"
    tlsextract XYZIN $pdbfile TLSOUT $tls_inp >/dev/null

    set size=`wc -l $tls_inp  | awk '{print $1}'`
    if( ! -e $tls_inp ||  $size <= 6 ) then
        /bin/rm -f $tls_inp
        echo "Error: TLS parameters were not extracted!"
        exit()
    endif
else  #tls.inp provided!
    set tls_inp = $2
endif


# do tlsanl to get a new pdb, memo below:
# BRESID [true | t | false | f]: if t, pdb contains B_partial; if f, B_full
# ISOOUT [FULL | RESI | TLSC]: RESI, B_partial; TLSC, TLS contribution

echo "Doing TLSANL to get full B factors($1)"
set tlspdb="${pdbfile}_tls"
set tlslog="${pdbfile}_tls.log"
tlsanl XYZIN  $pdbfile  TLSIN  $tls_inp  XYZOUT  $tlspdb <<end >$tlslog
BRESID %s
#NUMERIC  # a useful keyword, but only for new version (ccp4, 6.4)
ISOOUT %s 
end

""" %(bresid, isoout)
    pdb_tls = pdbfile + "_tls"
    log = pdb_tls + '.log'
    if os.path.exists(pdb_tls): os.remove(pdb_tls)
    if os.path.exists(log): os.remove(log) 

    script="tlsanl_script.csh"
    fo=open(script, "w")
    fo.write(csh_script)
    fo.close()
    command="chmod +x %s ; ./%s %s %s" % (script, script, pdbfile, tls_new_inp)
    os.system(command)

    return (pdb_tls ,log)


##########################################################
def delete_0_origin(pdbfile):
    '''delete the tls group with zero origins.This group is empty
    
    '''
    
    fp=open(pdbfile, 'r')
    pdb_new=pdbfile + '_000'
    fw=open(pdb_new, "w")
    i, idd, ntls, nn = 0, 0,0, []
    for x in fp:
        nn.append(x)
        if 'REMARK   3   TLS GROUP :' in x:
            ntls=x.split(':')[1].strip()
            
        elif 'ENDMDL' in x[:6] :
            print('This is a multi-model entry')
            #break
        
        elif 'REMARK   3    ORIGIN FOR THE GROUP (' in x and '):' in x :
            val=x.split(':')[1].split()
            if 'NULL' in val:
                print('Error: wrong origin ' + x )
                continue
            if float(x[39:48])==0.000 and  float(x[48:57])==0.000 :
                n=len(nn)
                for m in range(n): #move up (remove tls group)
                    k=n-m-1
                    #print(i, k, nn[k])
                    if 'REMARK   3   TLS GROUP :' in nn[k]: 
                        nn.pop(k)
                        break
                    nn.pop(k)

                    
                for x in fp: #move down
                    idd=idd+1
                    if 'S31:' in x: break
                    continue
            
        i=i+1
    for x in nn: fw.write(x)
    if(idd >3):
        t='Warning: tls group (%s) is deleted, since all parameters are zero.' %ntls
        perror(t)
    fp.close()
    fw.close()                
    return pdb_new 
    
##########################################################
def correct_tls(pdbfile):
    """ check PDB for possible corrections
    Correct TLS residue range if given as default (for TLSANL)
    """
    
    
    pdb_new=pdbfile + '_new'
    fw=open(pdb_new, "w")
    fr=open(pdbfile, "r")
    
    res_range=chain_res_range(pdbfile) #dic with chainID and residue range
    print('Polymer= ', res_range[1])
    print('Ligands= ', res_range[3])
    print('Waters = ', res_range[5])
    
    #print(res_range[2])
    #print(res_range[4])
    range1={}
    prog, remark_b, tls_b, anisou, ntls = '?'*5
    for line in fr:
        if('REMARK   3   PROGRAM     :' in line and 'REFMAC' in line):
            prog='REFMAC'
            
        elif "NUMBER OF TLS GROUPS  :" in line:
            ntls = get_value_after_id(line,":")
        elif "ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY" in line:
            tls_b = 'RESID'
        elif "ATOM RECORD CONTAINS SUM OF TLS AND RESIDUAL B " in line:
            tls_b = 'FULL'
        elif "REMARK   3  U VALUES      : RESIDUAL ONLY" in line:
            remark_b = 'RESID'
        elif "REMARK   3  U VALUES      : WITH TLS ADDED" in line:
            remark_b = 'FULL'
            
        elif "ANISOU" in line[0:6]:
            anisou = 'T'

        elif "REMARK   3    RESIDUE RANGE : " in line:
            tmp=line.split(':')[1].split()
            if (' NULL ' not in line and len(tmp) ==4 and tmp[0] == tmp[2]) :
                if is_digit(tmp[1])<=0 or is_digit(tmp[3])<=0 :
                   # print('Error: (%s; Residue-range: %s)' %(pdbfile, line[29:].strip()))
                    continue

                if tmp[0] not in range1.keys() and tmp[0] != ' ' : range1[tmp[0]]=[]
                range1[tmp[0]].append([int(eval(tmp[1])),int(eval(tmp[3]))])
                line=correct_tls_overlap(pdbfile, line, range1)
                if len(line) >0 : line=correct_tls_overlap(pdbfile, line, range1)
                if len(line)<5 :
                    print('Warning! (%s)skiping this line. \n' %(pdbfile))
                    #continue
            line1 = correct_residue_range(ntls, pdbfile, line, res_range)
            
            if len(line1) >10: fw.write(line1)
            continue
        
        elif line[:6] == 'ENDMDL' :
            print  ("Warning! The PDB file %s has multi-model. (1st model used)" % pdbfile)
            break       

        fw.write(line)
       
    fr.close()
    fw.close()
    if remark_b == 'FULL' or tls_b == 'FULL' or anisou == 'T' :
        print(pdbfile + "; Warning! PDB file with B_full (or with ANISOU record).\n")
#        sys.exit()
    return  pdb_new
        

##########################################################
def chain_res_range(pdbfile):
    '''use dic to contain chain-ID and residue range
    separate poly-peptide and hetatoms:
    '''

    pdb1, pdb2, pdb3 = separate_pdb(pdbfile)
    chain, chain_range = chain_res(pdb1,0) #poly-peptide
    chain1, chain_range1 = chain_res(pdb2,0) #ligands
    chain2, chain_range2 = chain_res(pdb3,1) #waters
    delete_file(pdb1, pdb2, pdb3)

    return chain, chain_range, chain1, chain_range1, chain2, chain_range2

##########################################################
def separate_pdb(pdbfile):
    '''separate PDB file by polymer pdb1 and ligands pdb2, water pdb3
    '''

    pdb1=pdbfile + 'tmp1'
    pdb2=pdbfile + 'tmp2'
    pdb3=pdbfile + 'tmp3'
    
    fw1, fw2, fw3 = open(pdb1,'w'), open(pdb2,'w'), open(pdb3,'w')
    
    fr, fr1 = [], open(pdbfile, "r").readlines()
    for x in fr1:
        if 'ANISOU' in x[:6] : continue
        if (atom_record(x) and ('HOH' in x[17:20] or 'DOD' in x[17:20])):
            fw3.write(x)
        else:
            fr.append(x)
 
    length = len(fr)    
    i, n = 0, 0
    for i in range(length):
        k=length-i-1
        if atom_record(fr[k]) :
            if residue(fr[k][17:20]) and 'HETA' not in fr[k][:4]:
                n=k 
                break
           
    if i> length-3 : n=length
    for i in range(length):
        if i<=n:
            fw1.write(fr[i])
        else:
            fw2.write(fr[i])

    fw1.close() 
    fw2.close() 
    fw3.close()
    return pdb1, pdb2, pdb3
##########################################################
def chain_res(pdbfile, id):
    '''use dic to contain chain-ID and residue range
    id=0: not include waters; id=1: include waters.
    '''
    
    fr=open(pdbfile, "r")
    
#    frn=[]
#    for ln in fr:
#        if('ATOM' in ln[:4] or 'HETATM' in ln[:6] ): frn.append(ln)
 #   
    chain={}
    chain_range={}
    n_old=-99999
    ch_old='?'
    for ln in fr:
        if('ATOM' in ln[:4] or 'HETATM' in ln[:6] ): 
            if id==0 and 'HOH' in ln[16:20] : continue
            ch=ln[20:22].strip()
            if ch not in chain.keys() and ch != ' ' : chain[ch]=[]
            nres=ln[22:26]
            n=int(nres)
            
            if ch in chain.keys():
                if(n == n_old and ch == ch_old) :
                    continue
                chain[ch].append(n)
            n_old=n
            ch_old=ch
            
            
    for key in chain:
        if(len(chain[key])>0): chain_range[key]= [min(chain[key]), max(chain[key])]

    fr.close()
    return chain, chain_range


##########################################################
def correct_tls_overlap(pdbfile, line, range):
    tmp=line.split(':')[1].split()
    ch, n1, n2 = tmp[0], int(eval(tmp[1])), int(eval(tmp[3]))  #from current line
    n=0
    ln=line
    #print(ch, range[ch])
    for x in reversed(range[ch]):
        n=n+1
        if n==1: continue
        k1, k2=x[0],x[1]
        #print(n1, n2, '|', k1, k2)
        if n2<n1 or k2<k1:
            print('Error: Residue range error [%d %d] ' %(n1, n2) , x)
            return ln
            
        elif (n1>=k1 and n2 <= k2) or (n1<k1 and n2 > k2) :
            print('Error: Total overlaps for [%d %d] ' %(n1, n2) , x)
            return ''
        
        elif (n1<=k2 and n1>=k1 and n2>k2 ):  #overlap current range with previous one
            t='%4s%6d%9s%6d' %(tmp[0],k2+1,tmp[0],n2)
            t1='Warning: residue range corrected.(%s)->(%s)' %(line[30:60].strip(),t)
            perror(t1)
            ln=line[:29] + t  + '   \n'
            break
        
        elif (n2>=k1 and n2<=k2 and n1<k1): #overlap current range with previous one
            t='%4s%6d%9s%6d' %(tmp[0],n1,tmp[0],k1-1)
            t1= 'Warning: residue range corrected.(%s)->(%s)' %(line[30:60].strip(),t)
            perror(t1)
            ln=line[:29] + t  + '   \n'
            break
        
        
    return ln

##########################################################
def correct_residue_range(ntls,pdbfile, line, res_range):
    '''Correct various errors in TLS residue range
    '''

    line1=line
    chain, chain_range, chain1, chain_range1, chain2, chain_range2 = res_range
    
    res=line[29:].split()
    nr = len(res)
    ch1, n1, ch2, n2 = (' ', 0, ' ', 0 )

    if 'NULL' in line or nr<2:
        print(pdbfile + '; Error: wrong residue range(%s).' %line[29:].strip())
    
    elif (nr==4): # most common cases
        ch1 = ch2 = res[0]
        if res[0] != res[2]:
            print(pdbfile + '; Error: wrong residue range(%s).\n' %line[29:].strip())

        else:
            line1 = correct_range(pdbfile, res[1],res[3], line, ch1, res_range)

    elif (nr==3):
        ch1 = ch2 = res[0]

        line1= correct_range(pdbfile, res[1],res[2], line, ch1, res_range)
        t='Warning: residue range correction(%s)->(%s).\n' %(line[29:55].strip(), line1[29:55].strip())
        perror(t)
        
    elif (nr==2 and is_digit(ntls)>0 and int(ntls)==1 and
          is_digit(res[0])>0 or is_digit(res[1])>0):
        if ('-' in res[0] and '-' in res[1]):
            line1 = ''
            if len(chain) >0 : line1=correct_default(line1,line, chain_range,  pdbfile)
            if len(chain1)>0 : line1=correct_default(line1,line, chain_range1, pdbfile)
            if len(chain2)>0 : line1=correct_default(line1,line, chain_range2, pdbfile)
            
        elif int(res[1])> int(res[0]) :
            ch1=list(chain.keys())[0]
            line1= correct_range(pdbfile, res[0],res[1], line, ch1, res_range)
            t='Warning: residue range correction(%s)->(%s).\n' %(line[29:55].strip(), line1[29:55].strip())
            perror(t)
            
        else :
            ch1=list(chain.keys())[0]
            line1= correct_range(pdbfile, res[0],res[1], line, ch1, res_range)

    return line1
   
##########################################################
def correct_default(line1, line, chain_range, pdbfile):
    for x,y in chain_range.items():
        rrange = '%4s%6d%9s%6d' %(x,y[0],x,y[1])
        tmp=line[:29] + rrange + '   \n'
        line1= tmp+line1
        t='Warning: residue range correction (%s)->(%s).\n' %(line[34:55].strip(), tmp[30:55].strip())
        perror(t)
    return line1
    
    
##########################################################
def correct_range(pdbfile, res1, res2, ln, ch, res_range):
    '''correct various residue range errors to pass the TLSANL
    '''

    chain,chain_range, chain1,chain_range1, chain2,chain_range2 = res_range

    line = ''

    if (is_digit(res1)==0 or is_digit(res2)==0)  :
        print(pdbfile + '; Error: wrong residue range(%s).' %ln[29:].strip())
        return ln
        
    elif (ch not in chain_range.keys() and ch not in chain_range1.keys() and
          ch not in chain_range2.keys()) :
        print(pdbfile + '; Error: Chain (%s) not in coordinate.' %ch)
        return ln
    
    n1, n2 = int(res1), int(res2) #tls residue number
    range1=[]
    if ch in chain_range.keys() :
        range1.append([chain_range[ch][0], chain_range[ch][1]])
    if ch in chain_range1.keys():
        range1.append([chain_range1[ch][0], chain_range1[ch][1]])
    if ch in chain_range2.keys():
        range1.append([chain_range2[ch][0], chain_range2[ch][1]])
        
    np=0
    for x in  range1: #
        if (n1 < x[0] and n2 < x[0]) or (n1 > x[1] and n2 > x[1]) : np = np + 1

    if (ch in chain_range.keys() and
        (n1 in chain[ch] or n2 in chain[ch] or
         (n1<chain_range[ch][0] and n2>chain_range[ch][1]) and n2<=999)): #polymer
        
        line=check_res_range(pdbfile, ln,chain, ch, n1, n2, 0)
        
    elif ch in chain_range1.keys() and (n1 in chain1[ch] or n2 in chain1[ch]): #ligands
        line=check_res_range(pdbfile, ln,chain1, ch, n1, n2, 1)
        
    elif ch in chain_range2.keys() and (n1 in chain2[ch] or n2 in chain2[ch]): #waters
        line=check_res_range(pdbfile, ln,chain2, ch, n1, n2, 2)
        
    elif ch in chain_range.keys() and ((n1<0 and n2<0) or (n1 <1 and  n2 >999))  : # default (all in ch)
        poly1, poly2 = chain_range[ch][0], chain_range[ch][1]
        tmp1, tmp2, tmp3 = '','',''
        if ch in chain_range.keys(): #poly
            n1, n2, idd = chain_range[ch][0], chain_range[ch][1], 1
            tmp1=make_new_line(pdbfile, ln, ch, n1, n2, 1)
            
        if ch in chain_range1.keys(): #lig
            n1, n2, idd = chain_range1[ch][0], chain_range1[ch][1], 1

            if n1<poly1 and n2>poly2: #lig residue range span polymer
                tmp2=make_two_newline(pdbfile, chain1[ch], ln, ch, poly1, poly2, 1)
            else:
                tmp2=make_new_line(pdbfile, ln, ch, n1, n2, 1)
            
        if ch in chain_range2.keys(): #water
            n1, n2, idd = chain_range2[ch][0], chain_range2[ch][1], 1
            
            if n1<poly1 and n2>poly2: #lig residue range span polymer
                tmp3=make_two_newline(pdbfile, chain2[ch], ln, ch, poly1, poly2, 1)
            else:
                tmp3=make_new_line(pdbfile, ln, ch, n1, n2, 1)
                
            
        line=tmp1 + tmp2 + tmp3
        
    elif (len(range1) == np and np>0 )  : # 
        line=' '
        t='Warning: Residue range [%s %d-%d] not exist in coordinate. This range removed.' %(ch,n1,n2)
        perror(t)
        
    if len(line)==0: line = ln

    return line

##########################################################
def check_res_range(pdbfile,ln, chain,  ch, n1, n2, idd):
    line=ln    
    if (n1<=n2) :
        idd, nd = 0, n2-n1
        if(n1 not in chain[ch]): n1=new_number(chain[ch],n1, n2, 1)
        if(n2 not in chain[ch]): n2=new_number(chain[ch],n1, n2,-1)

        if (nd != n2-n1): idd=1
        line = make_new_line(pdbfile, ln, ch, n1, n2, idd)
          
    elif (n1>n2) :
        line = make_new_line(pdbfile, ln, ch, n2, n1, 1)
    
            
    return line

##########################################################
def make_two_newline(pdbfile, resn_in, ln, ch, n1, n2, idd):
    
    resn=sorted(resn_in)
    #print(resn, ch, n1, n2, idd)
    k2, nn= 0, len(resn)-1
    for n in range( nn ):
        if resn[n+1] > n1 :
            k2=n
            break
        
    m1, m2, m3, m4 = resn[0], resn[k2],  resn[k2+1], resn[nn]
    #print(m1, m2, m3, m4, ';', resn, ch, n1, n2, idd)
    
    tmp1=make_new_line(pdbfile, ln, ch, m1, m2, idd)
    tmp2=make_new_line(pdbfile, ln, ch, m3, m4, idd)

    return tmp1+tmp2
        
    
##########################################################
def make_new_line(pdbfile, ln, ch, n1, n2, idd):
    ''' return a new line with correct residue ranges.
    id==1, corrected. idd==0, original
    '''
    
    tmp='%4s%6d%9s%6d' %(ch,n1,ch,n2)
    if idd>0:
        t='Warning: residue range correction(%s)->(%s).\n' %(ln[30:60].strip(),tmp.strip())
        perror(t)
    line1=ln[:29] + tmp  + '   \n'

    return line1
    
##########################################################
def new_number(chain, n1, n2, idd):

    if idd == 1: #increase the number
        n=n1
        for i in range(n2-n1+1):
            if n in chain or n>=n2 : break
            n=n1+i
        if(n-n1)>20:
            t='Warning: Residue range changed too much (%d).Watch out!' %(n-n1)
            perror(t)
       
    else:  #decrease the number
        n=n2
        for i in range(n2-n1+1):
            if n in chain or n<n1 : break
            n=n2-i
        
        if(n2-n)>20:
            t='Warning: Residue range changed too much (%d).Watch out!' %(n2-n)
            perror(t)
        
    return n

##########################################################
