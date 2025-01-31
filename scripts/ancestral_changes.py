import baltic as bt
import pandas as pd

treeFile = "annotated_tree.nexus"
outFile = 'annottated_tree_events.csv'
myTree=bt.loadNewick(treeFile, absoluteTime=False)
myTree.setAbsoluteTime(2023.923288) # need to set this to time of last sampled tip which for our data set is new date

myTree.traverse_tree() ## required to set heights
myTree.treeStats() ## report stats about tree

changes = 0
times = []
origins = []
destinations = []
for k in myTree.Objects: ## iterate over a flat list of branches

    "Assign node UNKNOWN region if not give"
    if 'region' in k.traits: # change to country or new_region if needed
        region = k.traits['region'] # change to country or new_region if needed
    else:
        region = 'UNKNOWN'
        k.traits['region'] = region # change to country or new_region if needed

    "Find parent region if given"
    if k.parent.traits:
        parent_region = k.parent.traits['region'] # change to country or new_region if needed
    else:
        parent_region = 'UNKNOWN'

    if region != parent_region:
        changes += 1
        times = times + [k.absoluteTime]
        origins = origins + [parent_region]
        destinations = destinations + [region]


print("Total number of state changes: " + str(changes))

df = pd.DataFrame({'EventTime' : times})
df['Origin'] = origins
df['Destination'] = destinations

df.to_csv(outFile)  
