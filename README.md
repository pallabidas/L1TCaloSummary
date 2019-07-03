# L1TCaloSummary

The package L1Trigger/L1TCaloSummary is prepared for monitoring the CMS Layer-1 Calorimeter Trigger.

It is a playpen for various tests.

Installation:

```bash
cmsrel CMSSW_10_6_0_pre4
cd CMSSW_10_6_0_pre4/src
cmsenv
git cms-init
git cms-merge-topic isobelojalvo:run3-dev-$CMSSW_VERSION
cd L1Trigger
git clone git@github.com:isobelojalvo/L1TCaloSummary.git
cd ..
scram b -j 8
```

Unit test for this directory can be run using:

```bash
pushd $CMSSW_BASE/src/L1Trigger/L1TCaloSummary/test;scram b runtests;popd
```

The CMSSW producer can be excercised using:

```bash
pushd $CMSSW_BASE/src/L1Trigger/L1TCaloSummary/test
cmsRun testL1TCaloSummary.py runNumber=260627 dataStream=/JetHT/Run2015D-v1/RAW
popd
```

Take a look at the resulting file /data/$USER/l1tCaloSummary-<runNumber>.root using root
