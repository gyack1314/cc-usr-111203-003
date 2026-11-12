"""显示已有原始摘要与配对检验报告，来源文件保持不动。"""
import argparse
import json
from pathlib import Path
from .io import load
from .domain import summarize
from .pairing import StructuralError
from .verify import run_verification


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['inspect', 'verify'])
    parser.add_argument('file')
    parser.add_argument('--out', help='检验报告输出路径（仅 verify），默认写到当前目录')
    a = parser.parse_args()
    if a.command == 'inspect':
        v, digest = load(a.file)
        print(json.dumps({'source_sha256': digest, 'summary': summarize(v)},
                         ensure_ascii=False, indent=2, allow_nan=False))
        return
    source = Path(a.file)
    try:
        report = run_verification(source)
    except StructuralError as e:
        raise SystemExit(f'资料违反配对契约，检验中止: {e}') from None
    if a.out:
        out = Path(a.out)
    else:
        out = Path.cwd() / f"verification-{source.stem}-{report['source']['sha256'][:8]}.json"
    if out.resolve() == source.resolve():
        raise SystemExit('输出路径不能与输入资料相同')
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False),
                   encoding='utf8')
    print(json.dumps({'written': str(out), 'source_sha256': report['source']['sha256'],
                      'pairing_counts': report['pairing']['counts']},
                     ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
