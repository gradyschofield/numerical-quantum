//
// Created by Grady Schofield on 5/12/21.
//

#ifndef ATOM_MATRIXELEMENT_H
#define ATOM_MATRIXELEMENT_H


class MatrixElement {
    int row;
    int column;

public:

    MatrixElement() = delete;

    MatrixElement(MatrixElement && e)
            : row(e.getRow()), column(e.getColumn())
    {
    }

    MatrixElement(MatrixElement const & e)
            : row(e.getRow()), column(e.getColumn())
    {
    }

    MatrixElement(int row, int column)
            : row(row), column(column)
    {
    }

    int getRow() const {
        return row;
    }

    int getColumn() const {
        return column;
    }

    bool operator==(MatrixElement const & e) const {
        return row == e.getRow() && column == e.getColumn();
    }
};

namespace std {
    template<> struct hash<MatrixElement>
    {
        size_t operator()(MatrixElement const & elem) const noexcept {
            return elem.getRow() ^ (elem.getColumn() << 1);
        }
    };
}


#endif //ATOM_MATRIXELEMENT_H
