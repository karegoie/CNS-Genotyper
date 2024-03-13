# CNS-Genotyper

## Cas9 Mutant Genotyper with Next Generation Sequencing!

Dependency:

- Bio = From reading fastq files to making alignments
- click = Supporting Command Line Interface
- xlsxwriter = Writing xlsx-shaped log file

## Data folder of InDel_Type_Counter:

To run the program, you must put the NGS data in here:

only file with name '.fastq' and '.fastq.gz' will be read by the program,

and other files like this will not be used as a data.

### main.py -x R2 -x undetermined 

To ignore some files after getting all data from the NGS raw result,

you can delete some of it manually, or using the config in the main.py code.

the config [-x R2] will ignore files with name 'R2', which mainly means Read 2.

# 2023-Winter Lab Internship

## TODO

ULTIMATE

- Raw 받아서 그냥 진행하는 것...으로 변경!!!

- 여러 개를 한 번에 seq... < 의미가 없을 수도 있다는 점. PCR 도중에 하나가 너무 Dominant 해질 수도.

- 병렬 처리보다는 안전한 직렬, 그런 형태로.

- Guide RNA seq가 하나.

01/05

- err 기준, '특정 indel 구간' 에 대해서...

- align 하고, err 값 찾기가... 아예 거기서 할까 싶기도 하고.

- 일단 둘 중에 높은 걸로 주기? 그건 좀 그렇다.

01/04

- Debug

01/03

- 49.txt homo 아닐 가능성. (오염이 그렇게 잘 되는 편?)


12/28

- 주요 에러: 길이가 딱 150으로 맞는 건가? d5: 뒤에 반드시 i5? 규칙을 알아야.

- Pseudocode 작성

- CLI...

- 'naming' 처리에 주의! 특히 괜히 복잡한 작업들...

12/27

- Pseudocode 작성에 조금 더 힘쓸 것,

- CLI 형태로, Interface 가능하게 하기,

- 아예 '특정 경로에 있는 모든 파일'을 이용하도록 할까? one Guide_RNA for all?

- CSV 생성: 여러 개의 파일 한꺼번에 진행 시 가능하도록.

12/26

- Pseudocode 및 실제 code에 대해,

- 아예 세 줄을 받아서 하는 '비교 함수'를 따로 작성하는 편이.

## Schedule

매주 화요일 4시 저널

## Memo

Pseudocode 먼저 쓰고 만들면 꽤 편한 것 같다! 특히 변수명의 컨셉까지.