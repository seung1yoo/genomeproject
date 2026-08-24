# genomeproject

WGS 결과 파일 인벤토리와 샘플별 gVCF 기반 특정 마커 빈도를 계산하는 CLI입니다. 원본 genomics 파일은 복사하거나 수정하지 않습니다.

## 입력

파일 manifest는 탭으로 구분된 long-format 파일입니다.

```tsv
sample_id	data_type	path
S001	cram	/data/wgs/S001/S001.cram
S001	gvcf	/data/wgs/S001/S001.g.vcf.gz
S001	vcf	/data/wgs/S001/S001.vcf.gz
```

`data_type`은 `cnv`, `cram`, `dtc`, `gvcf`, `sv`, `vcf` 중 하나입니다. 마커 파일 역시 TSV이며 한 행에 하나의 ALT를 기록합니다.

```tsv
marker_id	chrom	pos	ref	alt
marker-1	chr1	10001	A	G
```

마커 좌표와 대립유전자는 GRCh38 기준이어야 합니다. 추적되는 기본 설정을 직접 수정하지 말고 `configs/default.yaml`을 `configs/local.yaml`로 복사한 뒤 FASTA 경로를 지정하십시오. 로컬 설정과 실제 manifest·마커 파일은 Git에서 제외됩니다.

## 설치와 실행

```bash
conda env create -f environment.yml
conda activate genomeproject
cp configs/default.yaml configs/local.yaml

genomeproject inventory build \
  --manifest data/manifests/sample_files.tsv \
  --config configs/local.yaml \
  --run-id inventory-20260824

genomeproject markers validate \
  --markers data/markers/targets.tsv \
  --config configs/local.yaml

genomeproject markers frequency \
  --manifest data/manifests/sample_files.tsv \
  --markers data/markers/targets.tsv \
  --config configs/local.yaml \
  --run-id markers-20260824
```

`--skip-normalization`은 이미 정규화된 마커를 사용하거나 `bcftools`가 없는 개발 테스트에서만 사용합니다.

## 결과 해석

- `raw_*`: 완전한 GT가 있으면 GQ/DP와 무관하게 포함한 통계
- `qc_*`: 기본 `GQ >= 20`, `DP >= 8`을 통과한 genotype 통계
- `AC`: 대상 ALT allele 수
- `AN`: 실제 호출된 allele 수. 샘플 수의 두 배로 고정하지 않습니다.
- `AF`: `AC / AN`
- `audit_by_chrom/`: 샘플·마커별 GT와 판독/제외 사유

gVCF에 record가 없거나 no-call이면 homozygous reference로 간주하지 않습니다. `<NON_REF>` reference block이 실제 마커 위치를 포함하고 GT가 호출된 경우에만 대상 ALT dosage 0으로 계산합니다.
