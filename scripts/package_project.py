"""Create a reviewable local GitHub upload bundle without raw downloads or temporary files."""
import hashlib
import json
import platform
import zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    files=[]
    for pattern in ['README.md','.gitignore','requirements*.txt','data/*.csv','docs/**/*.md',
                    'scripts/*.py','results/**/*.json','results/**/*.md','results/**/*.csv',
                    'figures/*.png','figures/*.svg','output/pdf/*.pdf']:
        files.extend(ROOT.glob(pattern))
    files=sorted(set(p for p in files if p.is_file() and p.name!='release_manifest.json'))
    manifest={'release':'exploratory-v1.0','date':'2026-11-20','python':platform.python_version(),
              'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'excluded':'Original downloaded input files, caches and temporary PDF renderings. See README for data retrieval.',
              'status':'Local bundle only; not published to GitHub.'}
    target=ROOT/'results/release_manifest.json'
    target.write_text(json.dumps(manifest,indent=2)+'\n')
    files.append(target)
    out=ROOT/'output/packages/retina-rod-cone-exploratory-v1.0.zip'
    out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as archive:
        for p in files:
            archive.write(p,Path('retina-rod-cone-exploratory')/p.relative_to(ROOT))
    with zipfile.ZipFile(out) as archive:
        assert archive.testzip() is None
        assert not any('/data/raw/' in n or '/tmp/' in n or '__pycache__' in n for n in archive.namelist())
        print(f'{len(archive.namelist())} project files; ZIP integrity verified; {out.stat().st_size:,} bytes')
    print(out)

if __name__=='__main__': main()
