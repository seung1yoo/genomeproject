# genomeproject

WGS 결과 파일의 인벤토리를 만들고, 샘플별 gVCF에서 지정한 마커의 대립유전자 빈도를 계산하는 CLI입니다. 원본 genomics 파일은 읽기만 하며 복사하거나 수정하지 않습니다.

## 1. 설치

```bash
conda env create -f environment.yml
conda activate genomeproject
```

설치 확인:

```bash
genomeproject --help
pytest -q
```

DuckDB와 pysam은 오래된 Linux 컴파일러에서도 소스 빌드가 발생하지 않도록 conda-forge 바이너리로 설치합니다. 환경 생성이 중간에 실패했다면 환경을 제거한 뒤 다시 생성합니다.

```bash
conda env remove -n genomeproject
conda env create -f environment.yml
```

이미 clone한 저장소를 업데이트한 경우에는 다음을 실행합니다.

```bash
git pull origin main
pip install -e .
pytest -q
```

## 2. GRCh38 설정

추적되는 기본 설정을 직접 수정하지 말고 로컬 설정을 만듭니다.

```bash
cp configs/default.yaml configs/local.yaml
```

`configs/local.yaml`의 FASTA 경로를 실제 GRCh38 FASTA로 변경합니다.

```yaml
reference:
  build: GRCh38
  fasta: /path/to/Homo_sapiens_assembly38.fasta

inventory:
  expected_data_types: [cnv, cram, dtc, gvcf, sv, vcf]

markers:
  min_gq: 20
  min_dp: 8
  autosomes_only: true
  normalize: true
  workers: 4

output:
  root: results
```

FASTA와 samtools 인덱스가 모두 있어야 합니다.

```bash
ls -lh /path/to/Homo_sapiens_assembly38.fasta \
  /path/to/Homo_sapiens_assembly38.fasta.fai
```

`.fai`가 없다면 생성합니다.

```bash
samtools faidx /path/to/Homo_sapiens_assembly38.fasta
```

로컬 설정, 실제 manifest, 마커 파일, 결과 파일은 Git에서 제외됩니다.

## 3. 샘플 manifest 작성

manifest는 탭으로 구분된 long-format TSV입니다. 경로는 절대경로를 권장합니다.

```tsv
sample_id	data_type	path
S001	cram	/data/wgs/S001.cram
S001	gvcf	/data/wgs/S001.hard-filtered.gvcf.gz
```

`data_type`은 `cnv`, `cram`, `dtc`, `gvcf`, `sv`, `vcf` 중 하나입니다. 같은 샘플의 파일에는 같은 `sample_id`를 사용합니다.

### IPMI 디렉터리 예시

다음 구조처럼 CRAM과 gVCF 파일명이 sample ID로 시작하는 경우:

```text
/bc_storage/BioResources/IPMI/cram/22121301000110.cram
/bc_storage/BioResources/IPMI/gvcf/22121301000110.hard-filtered.gvcf.gz
```

manifest를 자동 생성할 수 있습니다.

```bash
mkdir -p data/manifests

{
  printf 'sample_id\tdata_type\tpath\n'

  find /bc_storage/BioResources/IPMI/cram \
    -type f -name "*.cram" \
    | sort \
    | awk -F/ '{
        sample=$NF
        sub(/\.cram$/, "", sample)
        print sample "\tcram\t" $0
      }'

  find /bc_storage/BioResources/IPMI/gvcf \
    -type f -name "*.gvcf.gz" \
    | sort \
    | awk -F/ '{
        sample=$NF
        sub(/\.hard-filtered\.gvcf\.gz$/, "", sample)
        print sample "\tgvcf\t" $0
      }'
} > data/manifests/sample_files.tsv
```

생성 결과와 행 수를 확인합니다.

```bash
head data/manifests/sample_files.tsv
wc -l data/manifests/sample_files.tsv
```

CRAM 또는 gVCF 중 한쪽만 있는 샘플을 찾습니다.

```bash
awk -F'\t' '
  NR > 1 {
    seen[$1 SUBSEP $2] = 1
    samples[$1] = 1
  }
  END {
    for (sample in samples) {
      cram = ((sample SUBSEP "cram") in seen)
      gvcf = ((sample SUBSEP "gvcf") in seen)
      if (!cram || !gvcf)
        print sample, "cram=" cram, "gvcf=" gvcf
    }
  }
' data/manifests/sample_files.tsv | sort
```

한쪽 파일만 존재하는 샘플을 manifest에서 억지로 제거할 필요는 없습니다. 인벤토리 결과에 누락으로 명시됩니다.

## 4. 파일 인벤토리 생성

```bash
genomeproject inventory build \
  --manifest data/manifests/sample_files.tsv \
  --config configs/local.yaml \
  --run-id inventory-cram-gvcf-YYYYMMDD
```

요약을 출력합니다.

```bash
genomeproject inventory summary \
  --catalog results/inventory-cram-gvcf-YYYYMMDD/inventory/inventory.duckdb
```

주요 결과 파일:

- `file_inventory.csv`/`.parquet`: 파일별 존재 여부, 크기, 인덱스, 문제 상태
- `inventory_summary.csv`/`.parquet`: 데이터 유형별 파일·샘플·누락·인덱스 요약
- `sample_completeness.csv`/`.parquet`: 샘플별 데이터 유형 보유 여부
- `inventory.duckdb`: 위 테이블을 포함한 조회용 catalog
- `run_config.yaml`: 실행 당시 해석된 설정과 입력 경로

인덱스가 없는 CRAM을 조회하려면:

```bash
python - <<'PY'
import duckdb

db = "results/inventory-cram-gvcf-YYYYMMDD/inventory/inventory.duckdb"
with duckdb.connect(db, read_only=True) as con:
    rows = con.execute("""
        SELECT sample_id, path, index_path, status
        FROM file_inventory
        WHERE data_type = 'cram'
          AND index_exists = false
        ORDER BY sample_id
    """).fetchall()

for row in rows:
    print(*row, sep="\t")
PY
```

heredoc의 마지막 `PY`는 공백 없이 줄 맨 앞에 있어야 합니다.

CRAM 상태는 인덱스를 만들기 전에 확인합니다.

```bash
samtools quickcheck -v /path/to/sample.cram
```

`missing EOF block`이 나오면 파일이 잘렸거나 불완전하게 복사됐을 수 있습니다. 바로 인덱싱하지 말고 원본 대조 또는 전체 읽기 검사를 먼저 수행합니다. 정상 CRAM의 인덱스는 다음처럼 생성할 수 있습니다.

```bash
samtools index -c /path/to/sample.cram
```

## 5. 마커 파일 작성 및 검증

마커 파일도 탭으로 구분된 TSV이며, 한 행에 ALT 하나만 기록합니다.

```tsv
marker_id	chrom	pos	ref	alt
marker-1	chr1	10001	A	G
```

규칙:

- 좌표와 REF/ALT는 GRCh38 기준이어야 합니다.
- 현재 기본 설정은 상염색체 1–22만 허용합니다.
- `3`과 `chr3`은 FASTA contig에 맞게 해석됩니다.
- multi-allelic rsID는 분석할 ALT마다 별도 행이 필요합니다.
- rsID는 locus 식별자이므로 REF/ALT를 별도로 확인해야 합니다.
- `marker_id` 및 `chrom:pos:ref:alt` 조합은 중복될 수 없습니다.

### ADIPOQ 마커 3개 예시

아래 값은 GRCh38의 VCF식 1-based 좌표입니다.

```bash
mkdir -p data/markers

printf 'marker_id\tchrom\tpos\tref\talt\nrs17300539\tchr3\t186841671\tG\tA\nrs1501299\tchr3\t186853334\tG\tT\nrs266729\tchr3\t186841685\tC\tG\n' \
  > data/markers/targets.tsv
```

FASTA의 REF와 bcftools 정규화를 검증합니다.

```bash
genomeproject markers validate \
  --markers data/markers/targets.tsv \
  --config configs/local.yaml
```

정상 출력 예시:

```json
{"markers": 3, "build": "GRCh38"}
```

`--skip-normalization`은 이미 정규화된 마커를 사용하거나 bcftools가 없는 개발 테스트에서만 사용합니다.

## 6. gVCF 마커 빈도 계산

수천 개의 gVCF를 처리할 때는 SSH 연결 종료에 대비해 tmux 또는 screen 안에서 실행하는 것을 권장합니다. conda로 tmux를 설치하려면:

```bash
conda install -c conda-forge tmux
tmux new -s genome-markers
```

분석 실행:

```bash
time genomeproject markers frequency \
  --manifest data/manifests/sample_files.tsv \
  --markers data/markers/targets.tsv \
  --config configs/local.yaml \
  --run-id adiponectin-markers-YYYYMMDD
```

tmux 세션에서 빠져나올 때는 `Ctrl+B` 다음 `D`, 다시 접속할 때는 다음을 사용합니다.

```bash
tmux attach -t genome-markers
```

주요 결과:

```text
results/adiponectin-markers-YYYYMMDD/
├── run_config.yaml
└── marker_frequency/
    ├── marker_frequency.csv
    ├── marker_frequency.parquet
    └── audit_by_chrom/
        └── chrom=chr3/
            ├── observations.csv
            └── observations.parquet
```

빈도 결과를 확인합니다.

```bash
column -s, -t \
  results/adiponectin-markers-YYYYMMDD/marker_frequency/marker_frequency.csv
```

샘플·마커별 상태와 QC 제외 사유를 집계합니다.

```bash
python - <<'PY'
import duckdb

path = "results/adiponectin-markers-YYYYMMDD/marker_frequency/audit_by_chrom/chrom=chr3/observations.parquet"
with duckdb.connect() as con:
    rows = con.execute("""
        SELECT status, qc_exclusion, count(*) AS n
        FROM read_parquet(?)
        GROUP BY status, qc_exclusion
        ORDER BY n DESC
    """, [path]).fetchall()

for row in rows:
    print(*row, sep="\t")
PY
```

## 7. 결과 해석

- `raw_*`: 완전한 GT가 있으면 GQ/DP와 무관하게 포함한 통계
- `qc_*`: 기본 `GQ >= 20`, `DP >= 8`을 통과한 genotype 통계
- `called_n`, `call_rate`: 호출된 샘플 수와 전체 cohort 대비 비율
- `AC`: 대상 ALT allele 수
- `AN`: 실제 호출된 allele 수로, 샘플 수의 두 배로 고정하지 않음
- `AF`: `AC / AN`
- `carrier_n`: 대상 ALT를 하나 이상 보유한 샘플 수
- `hom_ref_n`, `het_n`, `hom_alt_n`: 이배체 genotype별 샘플 수
- `status_counts`: 판독 상태별 샘플 수
- `qc_exclusion_counts`: 낮은 GQ/DP 등 QC 제외 사유별 샘플 수
- `audit_by_chrom/`: 샘플·마커별 GT, GQ, DP, dosage 및 판독/제외 사유

gVCF에 record가 없거나 no-call이면 homozygous reference로 간주하지 않습니다. `<NON_REF>` reference block이 실제 마커 위치를 포함하고 GT가 호출된 경우에만 대상 ALT dosage 0으로 계산합니다.

## 8. 인덱스 경고 처리

다음 경고는 실행을 중단시키는 오류는 아니지만 무시해서는 안 됩니다.

```text
[W::hts_idx_load3] The index file is older than the data file: sample.gvcf.gz.tbi
```

gVCF가 실제로 수정됐거나 파일 복사로 수정 시간만 바뀐 경우 모두 발생할 수 있습니다. 실행은 완료할 수 있지만 audit에서 `QUERY_ERROR`, `NO_RECORD` 등의 상태를 확인해야 합니다.

수정 시간이 오래된 `.tbi`를 찾으려면:

```bash
python - <<'PY'
from pathlib import Path

root = Path("/path/to/gvcf")
stale = []
for gvcf in root.glob("*.gvcf.gz"):
    index = Path(str(gvcf) + ".tbi")
    if index.is_file() and index.stat().st_mtime < gvcf.stat().st_mtime:
        stale.append((gvcf, index))

for gvcf, index in stale:
    print(gvcf, index, sep="\t")
print("stale_index_count:", len(stale))
PY
```

인덱스가 실제 데이터보다 오래된 것으로 확인되면 운영 저장소 정책과 쓰기 권한을 확인한 뒤 재생성합니다. `-f`는 기존 인덱스를 덮어씁니다.

```bash
tabix -f -p vcf /path/to/sample.hard-filtered.gvcf.gz
```
