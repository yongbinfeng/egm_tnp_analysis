import copy


## whoever wrote this code originally should take a few lessons in basic python operations
## this is, by all means, terrible....
## i'm not proud that i didn't fix it, i'm just lazy

def createBins( bining, cut ):
    binCut = None
    nbin = 0

    nbin = 1
    #index = range(len(bining))
    #for ix in range(len(index)):
    #    index[ix] = -1
    index = [-1 for i in bining]
    listOfIndex = []    
    listOfIndex.append( index )

    print('this is listOfIndex', listOfIndex)
    ### first map nD bins in a single list
    for iiv,iv in enumerate(bining.keys()): ##marc range(len(bining)):
    #for iv in range(len(bining)):
        var = bining[iv]['var']
        #python3if not bining[iv].has_key('type') or not bining[iv].has_key('bins'):
        if not 'type' in bining[iv] or not 'bins' in bining[iv]:
            print('bining is not complete for var {v}'.format(v=var))
            return listOfIndex
        nb1D = 1
        if   bining[iv]['type'] == 'float' :
            nb1D = len(bining[iv]['bins'])-1
        elif bining[iv]['type'] == 'int' :
            nb1D = len(bining[iv]['bins'])
        nbin = nbin * nb1D

        print('this is nb1D', nb1D)
        listOfIndexInit = copy.deepcopy(listOfIndex)
        for ib_v in range(nb1D):
            if ib_v == 0 :
                for ib in range(len(listOfIndex)):
                    listOfIndex[ib][iiv] = ib_v            
            else: 
                for ib in range(len(listOfIndexInit)):
                    listOfIndexInit[ib][iiv] = ib_v
                           
                listOfIndex.extend(copy.deepcopy(listOfIndexInit))

    listOfBins = []
    #ibin = 0
    nbins = len(listOfIndex)
    print('this is listOfIndex', listOfIndex)
    for ibin,ix in enumerate(listOfIndex):
        ### make bin definition
        binCut   = None
        binName  = 'bin%02d'%ibin
        if nbins > 100   :  binName  = 'bin%03d'%ibin
        if nbins > 1000  :  binName  = 'bin%04d'%ibin
        if nbins > 10000 :  binName  = 'bin%d'%ibin

        binTitle = ''
        binVars = {}
        if not cut is None:
            binCut = cut

        #print('binname', binName)

        for iv,v in enumerate(bining.keys()):
            var     = bining[v]['var']
            bins1D  = bining[v]['bins']
            varType = bining[v]['type']
            if varType == 'float' :
                if binCut is None: 
                    binCut   = '%s >= %f && %s < %f' % (var,bins1D[ix[iv]],var,bins1D[ix[iv]+1])
                    binTitle = '%1.3f < %s < %1.3f'  % (bins1D[ix[iv]],var,bins1D[ix[iv]+1])
                else:
                    binCut   = '%s && %s >= %f && %s < %f' % (binCut  ,var,bins1D[ix[iv]],var,bins1D[ix[iv]+1])
                    binTitle = '%s; %1.3f < %s < %1.3f'    % (binTitle,bins1D[ix[iv]],var,bins1D[ix[iv]+1])
                binName  = '%s_%s_%1.2fTo%1.2f'  % (binName ,var,bins1D[ix[iv]],bins1D[ix[iv]+1])
                binVars[var] = { 'min': bins1D[ix[iv]], 'max': bins1D[ix[iv]+1]}

                
            if varType == 'int' :
                if binCut is None: 
                    binCut   = '%s == %d' % (var,bins1D[ix[iv]])
                    binTitle = '%s = %d'  % (var,bins1D[ix[iv]])
                else:
                    binCut   = '%s && %s == %d' % (binCut,var,bins1D[ix[iv]])
                    binTitle = '%s; %s = %d'    % (binTitle,var,bins1D[ix[iv]])
                binName  = '%s_%sEq%d' % (binName ,var,bins1D[ix[iv]])
                binVars[var] = { 'min': bins1D[ix[iv]], 'max': bins1D[ix[iv]]}

            binName = binName.replace('-','m')
            binName = binName.replace('.','p')

        listOfBins.append({'cut' : binCut, 'title': binTitle, 'name' : binName, 'vars' : binVars })
 
    listOfVars = []
    for iv in  bining.keys(): ## marc range(len(bining)):
        listOfVars.append(bining[iv]['var'])
        
    binDefinition = {
        'vars' : listOfVars,
        'bins' : listOfBins
        } 
    return binDefinition

def tuneCuts( bindef, cuts ) :
    if cuts is None:
        return
    
    for ibin in cuts.keys():
        cut0 = bindef['bins'][ibin]['cut']
        cut1 = cuts[ibin]
        bindef['bins'][ibin]['cut'] = '%s && %s ' % (cut0,cut1)
    

