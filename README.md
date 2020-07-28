# Short recipe for CMGTools 

For the general recipe, [follow these instructions](https://twiki.cern.ch/twiki/bin/view/CMS/CMGToolsReleasesExperimental).

--------------

#### Set up CMSSW and the base git

```
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src
cmsenv
git cms-init
```

#### Add the central cmg-cmssw repository to get the Heppy 80X branch

```
git remote add cmg-central https://github.com/CERN-PH-CMG/cmg-cmssw.git -f -t heppy_104X_dev
```

#### Configure the sparse checkout, and get the base heppy packages

```
cp /home/ucl/cp3/dcransha/ttHdiff/sparse-checkout_104X_heppy .git/info/sparse-checkout
git checkout -b heppy_104X_dev cmg-central/heppy_104X_dev
```
N.B. If you don't have access to /home/ucl/cp3/dcransha/ttHdiff/sparse-checkout_104X_heppy, you can instead add the following to .git/info/sparse-checkout

```
/.gitignore
/.clang-tidy
/.clang-format
/PhysicsTools/Heppy/
/PhysicsTools/HeppyCore/
/PhysicsTools/JetMCAlgos/
/RecoTauTag/RecoTau/
/EgammaAnalysis/ElectronTools/
/RecoEgamma/EgammaTools/
/RecoEgamma/ElectronIdentification/
/RecoEgamma/PhotonIdentification/
```

#### get NanoAOD-tools (they are now used for processing inputs in NanoAOD format)
```
git clone https://github.com/sscruz/nanoAOD-tools.git -b postproc_dev PhysicsTools/NanoAODTools
```

#### Now get the CMGTools subsystem from the cmgtools-lite repository

```
git clone https://github.com/djcranshaw/cmgtools-lite.git -b 104X_dev_nano_tools CMGTools
cd CMGTools
```

#### Add your fork, and push the latest branch to it

```
git remote add origin  https://github.com/YOUR_GITHUB_REPOSITORY/cmgtools-lite.git
git push -u origin 104X_dev_nano_dev
```

#### Compile

```
cd $CMSSW_BASE/src
scram b -j 8
```
