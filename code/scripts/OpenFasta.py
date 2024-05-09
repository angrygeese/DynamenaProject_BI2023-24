import os

from _io import TextIOWrapper
from dataclasses import dataclass
from typing import Union


@dataclass
class FastaRecord:
    rec_id: str
    rec_desc: str
    rec_seq: str

    def __bool__(self):
        if self.rec_id or self.rec_desc:
            return True
        return False

    def __repr__(self):
        return f"ID: {self.rec_id}\nDescription: {self.rec_desc}\nSequence: {self.rec_seq}\n"


class OpenFasta:

    # nuc_alphabet = set('ATUGCNatugcn')
    # prot_alphabet = set('ACDEFGHIKLMNOPQRSTUVWY')

    def __init__(self, file: Union[str, os.PathLike], mode: str = "r"):
        self.file = file
        self.mode = mode
        self.handler = None
        self.current = None

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        if isinstance(file, (str, os.PathLike)):
            self._file = file
        else:
            raise ValueError(f"Expected str or os.PathLike object, not {type(file)}")

    def __enter__(self):
        self.file_open()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.file_close()

    def __iter__(self):
        self.file_open()
        return self

    def __next__(self):
        record = FastaRecord("", "", "")  # `record` хранит `FastaRecord`, полученный в текущей итерации `__next__`
        # логика такая: если текущая строка существует, входим в `while` и итерируемся по `self.handler`
        if self.current is None:
            raise TypeError("I/O operation not possible, open file first.")
        while isinstance(self.current, str):
            # если в строке есть ID или описание (строка начинается с '>'), то проверяем:
            if self.current.startswith(">"):
                # 1) мы дошли до неё в этом вызове `__next__` (`record` не пустой) - завершаем формирование `FastaRecord` и возвращаем его.
                if record:
                    record.rec_seq = "".join(record.rec_seq)
                    return record
                # 2) мы дошли до неё в предыдущем вызове `__next__` (`record` пустой) - создаём новый `FastaRecord`
                record_header = self.current.split(None, 1)  # https://github.com/biopython/biopython/blob/a45b4d003b34ff48b832d39abcb24202c8e32c8b/Bio/SeqIO/FastaIO.py#L201C42-L201C46
                if any(record_header) and len(record_header) == 1:
                    record_header.append('')
                rec_id, rec_desc = record_header
                record = FastaRecord(rec_id.lstrip(">"), rec_desc, [])
            # если в строке есть нуклеотиды - пишем их в `rec_seq` текущей `FastaRecord`
            # elif set(self.current).issubset(OpenFasta.nuc_alphabet) or set(self.current).issubset(OpenFasta.prot_alphabet):
            else:
                record.rec_seq.append(self.current)
            try:
                self.current = next(self.handler).strip()
            # если в ходе итерации возникает `StopIteration` - перезаписываем текущую строку как `-1`, завершаем формирование `FastaRecord` и возвращаем его.
            except StopIteration:
                self.current = -1
                record.rec_seq = "".join(record.rec_seq)
                return record
        # если текущая строка `-1` - то в цикле выше мы дошли до конца файла, поэтому при дальнейших попытках вызвать `__next__` вызываем `StopIteration`
        if self.current == -1:
            raise StopIteration

    def file_open(self):
        if not isinstance(self.handler, TextIOWrapper) or self.handler.closed:
            self.handler = open(self.file, mode=self.mode)
            self.current = next(self.handler).strip()

    def file_close(self):
        if isinstance(self.handler, TextIOWrapper) and not self.handler.closed:
            self.handler.close()

    def read_record(self):
        return next(self)

    def read_records(self):
        records = []
        next_record_exists = True
        while next_record_exists:
            try:
                next_record = next(self)
            except StopIteration:
                next_record_exists = False
            else:
                records.append(next_record)
        return records
