/***************************************************************************\
|* Function Parser for C++ v3.2.1                                          *|
|*-------------------------------------------------------------------------*|
|* Function optimizer                                                      *|
|*-------------------------------------------------------------------------*|
|* Copyright: Joel Yliluoma                                                *|
\***************************************************************************/

/* NOTE:
   This is a concatenation of all the header and source files of the
   original optimizer source code. All the code has been concatenated
   into this single file for convenience of usage (in other words, to
   simply use the optimizer, it's enough to add this file to the project
   rather than a multitude of files which the original optimizer source
   code is composed of).

   Thus this file exists for the usage of the Function parser library
   only, and is not suitable for developing it further. If you want to
   develop the library further, you should download the development
   version of the library, which has all the original source files.
 */

#include "io/fparser/fpconfig.h"
#ifdef FP_SUPPORT_OPTIMIZER

#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH

#include <stdio.h>

#ifdef _MSC_VER

typedef unsigned long long fphash_value_t;
#define FPHASH_CONST(x) x##ULL

#else

#include <stdint.h>
typedef uint_fast64_t fphash_value_t;
#define FPHASH_CONST(x) x##ULL

#endif

namespace FUNCTIONPARSERTYPES
{
    struct fphash_t
    {
        fphash_value_t hash1, hash2;

        bool operator==(const fphash_t& rhs) const
        { return hash1 == rhs.hash1 && hash2 == rhs.hash2; }

        bool operator!=(const fphash_t& rhs) const
        { return hash1 != rhs.hash1 || hash2 != rhs.hash2; }

        bool operator<(const fphash_t& rhs) const
        { return hash1 != rhs.hash1 ? hash1 < rhs.hash1 : hash2 < rhs.hash2; }
    };
}

#endif
#include <vector>
#include <utility>

#include "io/fparser/fpconfig.h"
#include "io/fparser/fparser.h"


#ifdef FP_EPSILON
 #define NEGATIVE_MAXIMUM (-FP_EPSILON)
#else
 #define NEGATIVE_MAXIMUM (-1e-14)
#endif

namespace FPoptimizer_Grammar
{
    struct Grammar;
}

namespace FPoptimizer_ByteCode
{
    class ByteCodeSynth;
}

namespace FPoptimizer_CodeTree
{
    class CodeTreeParserData;
    class CodeTree;

    class CodeTreeP
    {
    public:
        CodeTreeP()                   : p(0)   { }
        CodeTreeP(CodeTree*        b) : p(b)   { Birth(); }
        CodeTreeP(const CodeTreeP& b) : p(&*b) { Birth(); }

        inline CodeTree& operator* () const { return *p; }
        inline CodeTree* operator->() const { return p; }

        CodeTreeP& operator= (CodeTree*        b) { Set(b); return *this; }
        CodeTreeP& operator= (const CodeTreeP& b) { Set(&*b); return *this; }

        ~CodeTreeP() { Forget(); }

    private:
        inline static void Have(CodeTree* p2);
        inline void Forget();
        inline void Birth();
        inline void Set(CodeTree* p2);
    private:
        CodeTree* p;
    };

    class CodeTree
    {
        friend class CodeTreeParserData;
        friend class CodeTreeP;

        int RefCount;

    public:
        /* Describing the codetree node */
        unsigned Opcode;
        union
        {
            double   Value;   // In case of cImmed: value of the immed
            unsigned Var;     // In case of cVar:   variable number
            unsigned Funcno;  // In case of cFCall or cPCall
        };
        struct Param
        {
            CodeTreeP param; // param node
            bool      sign;  // true = negated or inverted

            Param()                           : param(),  sign()  {}
            Param(CodeTree*        p, bool s) : param(p), sign(s) {}
            Param(const CodeTreeP& p, bool s) : param(p), sign(s) {}

            bool IsIdenticalTo(const Param& b) const
            {
                return sign == b.sign && param->IsIdenticalTo(*b.param);
            }
        };

        // Parameters for the function
        //  These use the sign:
        //   For cAdd: operands to add together (0 to n)
        //             sign indicates that the value is negated before adding (0-x)
        //   For cMul: operands to multiply together (0 to n)
        //             sign indicates that the value is inverted before multiplying (1/x)
        //   For cAnd: operands to bitwise-and together (0 to n)
        //             sign indicates that the value is inverted before anding (!x)
        //   For cOr:  operands to bitwise-or together (0 to n)
        //             sign indicates that the value is inverted before orring (!x)
        //  These don't use the sign (sign is always false):
        //   For cMin: operands to select the minimum of
        //   For cMax: operands to select the maximum of
        //   For cImmed, not used
        //   For cVar,   not used
        //   For cIf:  operand 1 = condition, operand 2 = yes-branch, operand 3 = no-branch
        //   For anything else: the parameters required by the operation/function
        std::vector<Param> Params;

        /* Internal operation */
        FUNCTIONPARSERTYPES::fphash_t      Hash;
        size_t        Depth;
        CodeTree*     Parent;
        const FPoptimizer_Grammar::Grammar* OptimizedUsing;
    public:
        CodeTree();
        ~CodeTree();

        explicit CodeTree(double v); // produce an immed

        /* Generates a CodeTree from the given bytecode */
        static CodeTreeP GenerateFrom(
            const std::vector<unsigned>& byteCode,
            const std::vector<double>& immed,
            const FunctionParser::Data& data);

        void SynthesizeByteCode(
            std::vector<unsigned>& byteCode,
            std::vector<double>&   immed,
            size_t& stacktop_max);
        void SynthesizeByteCode(FPoptimizer_ByteCode::ByteCodeSynth& synth);

        /* Regenerates the hash.
         * child_triggered=false: Recurse to children
         * child_triggered=true:  Recurse to parents
         */
        void Rehash(bool child_triggered);
        void Recalculate_Hash_NoRecursion();

        void Sort();
        void Sort_Recursive();

        void SetParams(const std::vector<Param>& RefParams, bool do_clone = true);
        void AddParam(const Param& param);
        void DelParam(size_t index);

        /* Clones the tree. (For parameter duplication) */
        CodeTree* Clone();

        bool    IsImmed() const;
        double GetImmed() const { return Value; }
        bool    IsLongIntegerImmed() const { return IsImmed() && GetImmed() == (double)GetLongIntegerImmed(); }
        long   GetLongIntegerImmed() const { return (long)GetImmed(); }
        bool      IsVar() const;
        unsigned GetVar() const { return Var; }

        void NegateImmed() { if(IsImmed()) Value = -Value;       }
        void InvertImmed() { if(IsImmed()) Value = 1.0 / Value;  }
        void NotTheImmed() { if(IsImmed()) Value = Value == 0.0; }

        struct MinMaxTree
        {
            double min,max;
            bool has_min, has_max;
            MinMaxTree() : min(),max(),has_min(false),has_max(false) { }
            MinMaxTree(double mi,double ma): min(mi),max(ma),has_min(true),has_max(true) { }
            MinMaxTree(bool,double ma): min(),max(ma),has_min(false),has_max(true) { }
            MinMaxTree(double mi,bool): min(mi),max(),has_min(true),has_max(false) { }
        };
        /* This function calculates the minimum and maximum values
         * of the tree's result. If an estimate cannot be made,
         * has_min/has_max are indicated as false.
         */
        MinMaxTree CalculateResultBoundaries_do() const;
        MinMaxTree CalculateResultBoundaries() const;

        enum TriTruthValue { IsAlways, IsNever, Unknown };
        TriTruthValue GetEvennessInfo() const;

        bool IsAlwaysSigned(bool positive) const;
        bool IsAlwaysParity(bool odd) const
            { return GetEvennessInfo() == (odd?IsNever:IsAlways); }
        bool IsAlwaysInteger() const;

        void ConstantFolding();

        bool IsIdenticalTo(const CodeTree& b) const;

    private:
        CodeTree(const CodeTree&);
        CodeTree& operator=(const CodeTree&);
    };

    inline void CodeTreeP::Forget()
    {
        if(!p) return;
        p->RefCount -= 1;
        if(!p->RefCount) delete p;
        //assert(p->RefCount >= 0);
    }
    inline void CodeTreeP::Have(CodeTree* p2)
    {
        if(p2) ++(p2->RefCount);
    }
    inline void CodeTreeP::Birth()
    {
        Have(p);
    }
    inline void CodeTreeP::Set(CodeTree* p2)
    {
        Have(p2);
        Forget();
        p = p2;
    }
}
#define FPOPT_NAN_CONST (-1712345.25) /* Would use 0.0 / 0.0 here, but some compilers don't accept it. */

namespace FPoptimizer_CodeTree
{
    class CodeTree;
}

namespace FPoptimizer_Grammar
{
    typedef unsigned OpcodeType;

    enum TransformationType
    {
        None,    // default
        Negate,  // 0-x
        Invert,  // 1/x
        NotThe   // !x
    };

    enum ImmedConstraint_Value
    {
        ValueMask = 0x07,
        Value_AnyNum     = 0x0, // any value
        Value_EvenInt    = 0x1, // any even integer (0,2,4, etc)
        Value_OddInt     = 0x2, // any odd integer (1,3, etc)
        Value_IsInteger  = 0x3, // any integer-value (excludes e.g. 0.2)
        Value_NonInteger = 0x4  // any non-integer (excludes e.g. 1 or 5)
    };
    enum ImmedConstraint_Sign
    {
        SignMask  = 0x18,
        Sign_AnySign     = 0x00, // - or +
        Sign_Positive    = 0x08, // positive only
        Sign_Negative    = 0x10, // negative only
        Sign_NoIdea      = 0x18  // where sign cannot be guessed
    };
    enum ImmedConstraint_Oneness
    {
        OnenessMask   = 0x60,
        Oneness_Any      = 0x00,
        Oneness_One      = 0x20, // +1 or -1
        Oneness_NotOne   = 0x40  // anything but +1 or -1
    };

    enum SpecialOpcode
    {
        NumConstant = 0xFFFB, // Holds a particular value (syntax-time constant)
        ImmedHolder,          // Holds a particular immed
        NamedHolder,          // Holds a particular named param (of any kind)
        SubFunction,          // Holds an opcode and the params
        RestHolder            // Holds anything else
      //GroupFunction         // For parse-time functions
    };

    enum ParamMatchingType
    {
        PositionalParams, // this set of params in this order
        SelectedParams,   // this set of params in any order
        AnyParams         // these params are included
    };

    enum RuleType
    {
        ProduceNewTree, // replace self with the first (and only) from replaced_param
        ReplaceParams   // replace indicate params with replaced_params
    };

    enum SignBalanceType
    {
        BalanceDontCare,
        BalanceMoreNeg,
        BalanceMorePos,
        BalanceEqual
    };

    struct MatchResultType
    {
        //bool found:16;
        //bool has_more:16;
        // msvc doesn't like the above
        bool found;
        bool has_more;

        MatchResultType(bool f,bool m) : found(f),has_more(m) { }
    };
    static const MatchResultType
        NoMatch(false,false),       // No match, don't try to increment match_index
        TryMore(false,true),        // No match, but try to increment match_index
        FoundSomeMatch(true,true),  // Found match, but we may have more
        FoundLastMatch(true,false); // Found match, don't have more

    // For iterating through match candidates
    template<typename Payload>
    struct MatchPositionSpec
    {
        unsigned roundno;
        bool     done;
        Payload  data;
        MatchPositionSpec() : roundno(0), done(false), data() { }
    };

    /***/

#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif

    struct MatchedParams
    {
        ParamMatchingType type    : 4; // needs 2
        SignBalanceType   balance : 4; // needs 2
        // count,index to plist[]
        unsigned         count : 8;   // needs 2
        unsigned         index : 16;  // needs 10

        struct CodeTreeMatch;

        MatchResultType
            Match(FPoptimizer_CodeTree::CodeTree& tree,
                  CodeTreeMatch& match,
                  unsigned long match_index,
                  bool recursion) const;

        void ReplaceParams(FPoptimizer_CodeTree::CodeTree& tree,
                           const MatchedParams& matcher, CodeTreeMatch& match) const;

        void ReplaceTree(FPoptimizer_CodeTree::CodeTree& tree,
                         const MatchedParams& matcher, CodeTreeMatch& match) const;

        void SynthesizeTree(
            FPoptimizer_CodeTree::CodeTree& tree,
            const MatchedParams& matcher,
            MatchedParams::CodeTreeMatch& match) const;
    } PACKED_GRAMMAR_ATTRIBUTE;

    struct ParamSpec
    {
        OpcodeType opcode : 16;
        bool     sign     : 1;
        TransformationType
           transformation  : 3;
        unsigned minrepeat : 3;
        bool     anyrepeat : 1;

        // For NumConstant:   index to clist[]
        // For ImmedHolder:   index is the slot, count is an ImmedConstraint
        // For RestHolder:    index is the slot
        // For NamedHolder:   index is the slot, count is an ImmedConstraint
        // For SubFunction:   index to flist[],  count is an ImmedConstraint
        // For anything else
        //  =  GroupFunction: index,count to plist[]
        unsigned count : 8;
        unsigned index : 16;

        MatchResultType Match(
            FPoptimizer_CodeTree::CodeTree& tree,
            MatchedParams::CodeTreeMatch& match,
            TransformationType transf,
            unsigned long match_index) const;

        bool GetConst(
            const MatchedParams::CodeTreeMatch& match,
            double& result) const;

        void SynthesizeTree(
            FPoptimizer_CodeTree::CodeTree& tree,
            const MatchedParams& matcher,
            MatchedParams::CodeTreeMatch& match) const;
    } PACKED_GRAMMAR_ATTRIBUTE;
    struct Function
    {
        OpcodeType opcode : 16;
        // index to mlist[]
        unsigned   index  : 16;

        MatchResultType
            Match(FPoptimizer_CodeTree::CodeTree& tree,
                  MatchedParams::CodeTreeMatch& match,
                  unsigned long match_index) const;
    } PACKED_GRAMMAR_ATTRIBUTE;
    struct Rule
    {
        unsigned  n_minimum_params : 8;
        RuleType  type             : 8;
        // index to mlist[]
        unsigned  repl_index       : 16;
        Function  func;

        bool ApplyTo(FPoptimizer_CodeTree::CodeTree& tree) const;
    } PACKED_GRAMMAR_ATTRIBUTE;
    struct Grammar
    {
        // count,index to rlist[]
        unsigned index : 16;
        unsigned count : 16;

        bool ApplyTo(FPoptimizer_CodeTree::CodeTree& tree,
                     bool recursion=false) const;
    } PACKED_GRAMMAR_ATTRIBUTE;

    extern const struct GrammarPack
    {
        const double*         clist;
        const ParamSpec*      plist;
        const MatchedParams*  mlist;
        const Function*       flist;
        const Rule*           rlist;
        Grammar               glist[4];
    } pack;
}
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define CONSTANT_E     2.71828182845904509080  // exp(1)
#define CONSTANT_EI    0.3678794411714423215955 // exp(-1)
#define CONSTANT_2E    7.3890560989306502272304 // exp(2)
#define CONSTANT_2EI   0.135335283236612691894 // exp(-2)
#define CONSTANT_PI    M_PI                    // atan2(0,-1)
#define CONSTANT_L10   2.30258509299404590109  // log(10)
#define CONSTANT_L2    0.69314718055994530943  // log(2)
#define CONSTANT_L10I  0.43429448190325176116  // 1/log(10)
#define CONSTANT_L2I   1.4426950408889634074   // 1/log(2)
#define CONSTANT_L10E  CONSTANT_L10I           // log10(e)
#define CONSTANT_L10EI CONSTANT_L10            // 1/log10(e)
#define CONSTANT_L10B  0.3010299956639811952137 // log10(2)
#define CONSTANT_L10BI 3.3219280948873623478703 // 1/log10(2)
#define CONSTANT_LB10  CONSTANT_L10BI          // log2(10)
#define CONSTANT_LB10I CONSTANT_L10B           // 1/log2(10)
#define CONSTANT_L2E   CONSTANT_L2I            // log2(e)
#define CONSTANT_L2EI  CONSTANT_L2             // 1/log2(e)
#define CONSTANT_DR    (180.0 / M_PI)          // 180/pi
#define CONSTANT_RD    (M_PI / 180.0)          // pi/180

#define CONSTANT_POS_INF     HUGE_VAL  // positive infinity, from math.h
#define CONSTANT_NEG_INF   (-HUGE_VAL) // negative infinity
#define CONSTANT_PIHALF (M_PI / 2)
/* crc32 */

#ifdef _MSC_VER

 typedef unsigned int crc32_t;

#else

 #include <stdint.h>
 typedef uint_least32_t crc32_t;

#endif

namespace crc32
{
    enum { startvalue = 0xFFFFFFFFUL, poly = 0xEDB88320UL };

    /* This code constructs the CRC32 table at compile-time,
     * avoiding the need for a huge explicitly written table of magical numbers. */
    template<crc32_t crc> // One byte of a CRC32 (eight bits):
    struct b8
    {
        enum { b1 = (crc & 1) ? (poly ^ (crc >> 1)) : (crc >> 1),
               b2 = (b1  & 1) ? (poly ^ (b1  >> 1)) : (b1  >> 1),
               b3 = (b2  & 1) ? (poly ^ (b2  >> 1)) : (b2  >> 1),
               b4 = (b3  & 1) ? (poly ^ (b3  >> 1)) : (b3  >> 1),
               b5 = (b4  & 1) ? (poly ^ (b4  >> 1)) : (b4  >> 1),
               b6 = (b5  & 1) ? (poly ^ (b5  >> 1)) : (b5  >> 1),
               b7 = (b6  & 1) ? (poly ^ (b6  >> 1)) : (b6  >> 1),
               res= (b7  & 1) ? (poly ^ (b7  >> 1)) : (b7  >> 1) };
    };
    inline crc32_t update(crc32_t crc, unsigned/* char */b) // __attribute__((pure))
    {
        // Four values of the table
        #define B4(n) b8<n>::res,b8<n+1>::res,b8<n+2>::res,b8<n+3>::res
        // Sixteen values of the table
        #define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
        // The whole table, index by steps of 16
        static const crc32_t table[256] =
        { R(0x00),R(0x10),R(0x20),R(0x30), R(0x40),R(0x50),R(0x60),R(0x70),
          R(0x80),R(0x90),R(0xA0),R(0xB0), R(0xC0),R(0xD0),R(0xE0),R(0xF0) };
        #undef R
        #undef B4
        return ((crc >> 8) /* & 0x00FFFFFF*/) ^ table[/*(unsigned char)*/(crc^b)&0xFF];
    }
    inline crc32_t calc_upd(crc32_t c, const unsigned char* buf, size_t size)
    {
        crc32_t value = c;
        for(size_t p=0; p<size; ++p) value = update(value, buf[p]);
        return value;
    }
    inline crc32_t calc(const unsigned char* buf, size_t size)
    {
        return calc_upd(startvalue, buf, size);
    }
}
#include <string>

const std::string FP_GetOpcodeName(unsigned opcode, bool pad=false);
#include <string>
#include <sstream>
#include <assert.h>

#include <iostream>

#include "io/fparser/fpconfig.h"
#include "io/fparser/fptypes.h"


using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

const std::string FP_GetOpcodeName(unsigned opcode, bool pad)
{
#if 1
    /* Symbolic meanings for the opcodes? */
    const char* p = 0;
    switch(OPCODE(opcode))
    {
        case cAbs: p = "cAbs"; break;
        case cAcos: p = "cAcos"; break;
        case cAcosh: p = "cAcosh"; break;
        case cAsin: p = "cAsin"; break;
        case cAsinh: p = "cAsinh"; break;
        case cAtan: p = "cAtan"; break;
        case cAtan2: p = "cAtan2"; break;
        case cAtanh: p = "cAtanh"; break;
        case cCeil: p = "cCeil"; break;
        case cCos: p = "cCos"; break;
        case cCosh: p = "cCosh"; break;
        case cCot: p = "cCot"; break;
        case cCsc: p = "cCsc"; break;
        case cEval: p = "cEval"; break;
        case cExp: p = "cExp"; break;
        case cExp2: p = "cExp2"; break;
        case cFloor: p = "cFloor"; break;
        case cIf: p = "cIf"; break;
        case cInt: p = "cInt"; break;
        case cLog: p = "cLog"; break;
        case cLog2: p = "cLog2"; break;
        case cLog10: p = "cLog10"; break;
        case cMax: p = "cMax"; break;
        case cMin: p = "cMin"; break;
        case cPow: p = "cPow"; break;
        case cSec: p = "cSec"; break;
        case cSin: p = "cSin"; break;
        case cSinh: p = "cSinh"; break;
        case cSqrt: p = "cSqrt"; break;
        case cTan: p = "cTan"; break;
        case cTanh: p = "cTanh"; break;
        case cImmed: p = "cImmed"; break;
        case cJump: p = "cJump"; break;
        case cNeg: p = "cNeg"; break;
        case cAdd: p = "cAdd"; break;
        case cSub: p = "cSub"; break;
        case cMul: p = "cMul"; break;
        case cDiv: p = "cDiv"; break;
        case cMod: p = "cMod"; break;
        case cEqual: p = "cEqual"; break;
        case cNEqual: p = "cNEqual"; break;
        case cLess: p = "cLess"; break;
        case cLessOrEq: p = "cLessOrEq"; break;
        case cGreater: p = "cGreater"; break;
        case cGreaterOrEq: p = "cGreaterOrEq"; break;
        case cNot: p = "cNot"; break;
        case cAnd: p = "cAnd"; break;
        case cOr: p = "cOr"; break;
        case cDeg: p = "cDeg"; break;
        case cRad: p = "cRad"; break;
        case cFCall: p = "cFCall"; break;
        case cPCall: p = "cPCall"; break;
#ifdef FP_SUPPORT_OPTIMIZER
        case cVar: p = "cVar"; break;
        case cDup: p = "cDup"; break;
        case cInv: p = "cInv"; break;
        case cFetch: p = "cFetch"; break;
        case cPopNMov: p = "cPopNMov"; break;
        case cSqr: p = "cSqr"; break;
        case cRDiv: p = "cRDiv"; break;
        case cRSub: p = "cRSub"; break;
        case cNotNot: p = "cNotNot"; break;
        case cRSqrt: p = "cRSqrt"; break;
#endif
        case cNop: p = "cNop"; break;
        case VarBegin: p = "VarBegin"; break;
    }
    switch( SpecialOpcode(opcode) )
    {
        case NumConstant:   p = "NumConstant"; break;
        case ImmedHolder:   p = "ImmedHolder"; break;
        case NamedHolder:   p = "NamedHolder"; break;
        case RestHolder:    p = "RestHolder"; break;
        case SubFunction:   p = "SubFunction"; break;
      //case GroupFunction: p = "GroupFunction"; break;
    }
    std::ostringstream tmp;
    //if(!p) std::cerr << "o=" << opcode << "\n";
    assert(p);
    tmp << p;
    if(pad) while(tmp.str().size() < 12) tmp << ' ';
    return tmp.str();
#else
    /* Just numeric meanings */
    std::ostringstream tmp;
    tmp << opcode;
    if(pad) while(tmp.str().size() < 5) tmp << ' ';
    return tmp.str();
#endif
}
#include <vector>
#include <utility>

#include "io/fparser/fpconfig.h"
#include "io/fparser/fparser.h"
#include "io/fparser/fptypes.h"


namespace FPoptimizer_ByteCode
{
    class ByteCodeSynth
    {
    public:
        ByteCodeSynth()
            : ByteCode(), Immed(), StackTop(0), StackMax(0)
        {
            /* estimate the initial requirements as such */
            ByteCode.reserve(64);
            Immed.reserve(8);
        }

        void Pull(std::vector<unsigned>& bc,
                  std::vector<double>&   imm,
                  size_t& StackTop_max)
        {
            ByteCode.swap(bc);
            Immed.swap(imm);
            StackTop_max = StackMax;
        }

        size_t GetByteCodeSize() const { return ByteCode.size(); }
        size_t GetStackTop()     const { return StackTop; }

        void PushVar(unsigned varno)
        {
            ByteCode.push_back(varno);
            SetStackTop(StackTop+1);
        }

        void PushImmed(double immed)
        {
            using namespace FUNCTIONPARSERTYPES;
            ByteCode.push_back(cImmed);
            Immed.push_back(immed);
            SetStackTop(StackTop+1);
        }

        void StackTopIs(FUNCTIONPARSERTYPES::fphash_t hash)
        {
            if(StackTop > 0)
            {
                StackHash[StackTop-1].first = true;
                StackHash[StackTop-1].second = hash;
            }
        }

        void AddOperation(unsigned opcode, unsigned eat_count, unsigned produce_count = 1)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop - eat_count);

            if(opcode == cMul && ByteCode.back() == cDup)
                ByteCode.back() = cSqr;
            else
                ByteCode.push_back(opcode);
            SetStackTop(StackTop + produce_count);
        }

        void DoPopNMov(size_t targetpos, size_t srcpos)
        {
            using namespace FUNCTIONPARSERTYPES;
            ByteCode.push_back(cPopNMov);
            ByteCode.push_back( (unsigned) targetpos);
            ByteCode.push_back( (unsigned) srcpos);

            SetStackTop(srcpos+1);
            StackHash[targetpos] = StackHash[srcpos];
            SetStackTop(targetpos+1);
        }

        void DoDup(size_t src_pos)
        {
            using namespace FUNCTIONPARSERTYPES;
            if(src_pos == StackTop-1)
            {
                ByteCode.push_back(cDup);
            }
            else
            {
                ByteCode.push_back(cFetch);
                ByteCode.push_back( (unsigned) src_pos);
            }
            SetStackTop(StackTop + 1);
            StackHash[StackTop-1] = StackHash[src_pos];
        }

        bool FindAndDup(FUNCTIONPARSERTYPES::fphash_t hash)
        {
            for(size_t a=StackHash.size(); a-->0; )
            {
                if(StackHash[a].first && StackHash[a].second == hash)
                {
                    DoDup(a);
                    return true;
                }
            }
            return false;
        }

        void SynthIfStep1(size_t& ofs)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // the If condition was popped.

            ofs = ByteCode.size();
            ByteCode.push_back(cIf);
            ByteCode.push_back(0); // code index
            ByteCode.push_back(0); // Immed index
        }
        void SynthIfStep2(size_t& ofs)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // ignore the pushed then-branch result.

            ByteCode[ofs+1] = unsigned( ByteCode.size()+2 );
            ByteCode[ofs+2] = unsigned( Immed.size()      );

            ofs = ByteCode.size();
            ByteCode.push_back(cJump);
            ByteCode.push_back(0); // code index
            ByteCode.push_back(0); // Immed index
        }
        void SynthIfStep3(size_t& ofs)
        {
            SetStackTop(StackTop-1); // ignore the pushed else-branch result.

            ByteCode[ofs+1] = unsigned( ByteCode.size()-1 );
            ByteCode[ofs+2] = unsigned( Immed.size()      );

            SetStackTop(StackTop+1); // one or the other was pushed.
        }

    private:
        void SetStackTop(size_t value)
        {
            StackTop = value;
            if(StackTop > StackMax) StackMax = StackTop;
            StackHash.resize(value);
        }

    private:
        std::vector<unsigned> ByteCode;
        std::vector<double>   Immed;

        std::vector<
            std::pair<bool/*known*/, FUNCTIONPARSERTYPES::fphash_t/*hash*/>
                   > StackHash;
        size_t StackTop;
        size_t StackMax;
    };

    struct SequenceOpCode;
    extern const SequenceOpCode AddSequence; /* Multiplication implemented with adds */
    extern const SequenceOpCode MulSequence; /* Exponentiation implemented with muls */

    /* Generate a sequence that multiplies or exponentifies the
     * last operand in the stack by the given constant integer
     * amount (positive or negative).
     */
    void AssembleSequence(
        long count,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth);
}

using namespace FUNCTIONPARSERTYPES;

namespace FPoptimizer_ByteCode
{
    const struct SequenceOpCode
    {
        double basevalue;
        unsigned op_flip;
        unsigned op_normal, op_normal_flip;
        unsigned op_inverse, op_inverse_flip;
    } AddSequence = {0.0, cNeg, cAdd, cAdd, cSub, cRSub },
      MulSequence = {1.0, cInv, cMul, cMul, cDiv, cRDiv };
}

using namespace FPoptimizer_ByteCode;

#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
#ifndef FP_GENERATING_POWI_TABLE
static const
#endif
signed char powi_table[POWI_TABLE_SIZE] =
{
      0,   1,   1,   1,   2,   1,   3,   1, /*   0 -   7 */
      4,   1,   5,   1,   6,   1,   7,   5, /*   8 -  15 */
      8,   1,   9,   1,  10,   1,  11,   1, /*  16 -  23 */
     12,   5,  13,   9,  14,   1,  15,   1, /*  24 -  31 */
     16,   1,  17,   1,  18,   1,  19,  13, /*  32 -  39 */
     20,   1,  21,   1,  22,   9,   1,   2, /*  40 -  47 */
     24,   1,  25,  17,  26,   1,  27,  11, /*  48 -  55 */
     28,  19,  29,   8,  30,   1,  31,  21, /*  56 -  63 */
     32,   1,  33,   1,  34,   1,  35,   1, /*  64 -  71 */
     36,   1,  37,  25,  38,   1,  39,   1, /*  72 -  79 */
     40,   9,  41,   1,  42,  17,   1,  29, /*  80 -  87 */
     44,   1,  45,   1,  46,  31,  47,  19, /*  88 -  95 */
     48,   1,  49,  33,  50,   1,  51,   1, /*  96 - 103 */
     52,  35,  53,   8,  54,   1,  55,  37, /* 104 - 111 */
     56,   1,  57,  16,  58,  13,  59,  17, /* 112 - 119 */
     60,   1,  61,  41,  62,  25,  63,   1, /* 120 - 127 */
     64,   1,  65,   1,  66,   1,  67,  45, /* 128 - 135 */
     68,   1,  69,   1,  70,   1,  71,   8, /* 136 - 143 */
     72,   1,  73,  49,  74,   1,  75,   1, /* 144 - 151 */
     76,  17,   1,  31,  78,   1,  79,  53, /* 152 - 159 */
     80,   1,  81,   1,  82,  33,   1,   2, /* 160 - 167 */
     84,   1,  85,  19,  86,   8,  87,  35, /* 168 - 175 */
     88,   1,  89,   1,  90,   1,  91,  61, /* 176 - 183 */
     92,  37,  93,  17,  94,  21,  95,   1, /* 184 - 191 */
     96,   1,  97,  65,  98,   1,  99,   1, /* 192 - 199 */
    100,  67, 101,   8, 102,  41, 103,  69, /* 200 - 207 */
    104,   1, 105,  16, 106,  71, 107,   1, /* 208 - 215 */
    108,   1, 109,  73, 110,  17, 111,   1, /* 216 - 223 */
    112,  45, 113,  32, 114,   1, 115,  33, /* 224 - 231 */
    116,   1, 117,   1, 118,   1, 119,   1, /* 232 - 239 */
    120,   1, 121,  81, 122,  49, 123,  19, /* 240 - 247 */
    124,   1, 125,   1, 126,   1, 127,  85  /* 248 - 255 */
}; /* as in gcc, but custom-optimized for stack calculation */
static const int POWI_CACHE_SIZE = 256;

#define FPO(x) /**/
//#define FPO(x) x


namespace
{
    class PowiCache
    {
    private:
        int cache[POWI_CACHE_SIZE];
        int cache_needed[POWI_CACHE_SIZE];

    public:
        PowiCache()
            : cache(), cache_needed() /* Assume we have no factors in the cache */
        {
            /* Decide which factors we would need multiple times.
             * Output:
             *   cache[]        = these factors were generated
             *   cache_needed[] = number of times these factors were desired
             */
           for (size_t i=0; i<POWI_CACHE_SIZE; i++)
              cache[i] = 0;
           cache[1] = 1; // We have this value already.
        }

        bool Plan_Add(long value, int count)
        {
            if(value >= POWI_CACHE_SIZE) return false;
            //FPO(fprintf(stderr, "%ld will be needed %d times more\n", count, need_count));
            cache_needed[value] += count;
            return cache[value] != 0;
        }

        void Plan_Has(long value)
        {
            if(value < POWI_CACHE_SIZE)
                cache[value] = 1; // This value has been generated
        }

        void Start(size_t value1_pos)
        {
            for(int n=2; n<POWI_CACHE_SIZE; ++n)
                cache[n] = -1; /* Stack location for each component */

            Remember(1, value1_pos);

            DumpContents();
        }

        int Find(long value) const
        {
            if(value < POWI_CACHE_SIZE)
            {
                if(cache[value] >= 0)
                {
                    // found from the cache
                    FPO(fprintf(stderr, "* I found %ld from cache (%u,%d)\n",
                        value, (unsigned)cache[value], cache_needed[value]));
                    return cache[value];
                }
            }
            return -1;
        }

        void Remember(long value, size_t stackpos)
        {
            if(value >= POWI_CACHE_SIZE) return;

            FPO(fprintf(stderr, "* Remembering that %ld can be found at %u (%d uses remain)\n",
                value, (unsigned)stackpos, cache_needed[value]));
            cache[value] = (int) stackpos;
        }

        void DumpContents() const
        {
            FPO(for(int a=1; a<POWI_CACHE_SIZE; ++a)
                if(cache[a] >= 0 || cache_needed[a] > 0)
                {
                    fprintf(stderr, "== cache: sp=%d, val=%d, needs=%d\n",
                        cache[a], a, cache_needed[a]);
                })
        }

        int UseGetNeeded(long value)
        {
            if(value >= 0 && value < POWI_CACHE_SIZE)
                return --cache_needed[value];
            return 0;
        }
    };


    size_t AssembleSequence_Subdivide(
        long count,
        PowiCache& cache,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth);

    void Subdivide_Combine(
        size_t apos, long aval,
        size_t bpos, long bval,
        PowiCache& cache,

        unsigned cumulation_opcode,
        unsigned cimulation_opcode_flip,

        ByteCodeSynth& synth);

    void PlanNtimesCache
        (long value,
         PowiCache& cache,
         int need_count,
         int recursioncount=0)
    {
        if(value < 1) return;

    #ifdef FP_GENERATING_POWI_TABLE
        if(recursioncount > 32) throw false;
    #endif

        if(cache.Plan_Add(value, need_count)) return;

        long half = 1;
        if(value < POWI_TABLE_SIZE)
            half = powi_table[value];
        else if(value & 1)
            half = value & ((1 << POWI_WINDOW_SIZE) - 1); // that is, value & 7
        else
            half = value / 2;

        long otherhalf = value-half;
        if(half > otherhalf || half<0) std::swap(half,otherhalf);

        FPO(fprintf(stderr, "value=%ld, half=%ld, otherhalf=%ld\n", value,half,otherhalf));

        if(half == otherhalf)
        {
            PlanNtimesCache(half,      cache, 2, recursioncount+1);
        }
        else
        {
            PlanNtimesCache(half,      cache, 1, recursioncount+1);
            PlanNtimesCache(otherhalf>0?otherhalf:-otherhalf,
                                       cache, 1, recursioncount+1);
        }

        cache.Plan_Has(value);
    }

    size_t AssembleSequence_Subdivide(
        long value,
        PowiCache& cache,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth)
    {
        int cachepos = cache.Find(value);
        if(cachepos >= 0)
        {
            // found from the cache
            return cachepos;
        }

        long half = 1;
        if(value < POWI_TABLE_SIZE)
            half = powi_table[value];
        else if(value & 1)
            half = value & ((1 << POWI_WINDOW_SIZE) - 1); // that is, value & 7
        else
            half = value / 2;
        long otherhalf = value-half;
        if(half > otherhalf || half<0) std::swap(half,otherhalf);

        FPO(fprintf(stderr, "* I want %ld, my plan is %ld + %ld\n", value, half, value-half));

        if(half == otherhalf)
        {
            size_t half_pos = AssembleSequence_Subdivide(half, cache, sequencing, synth);

            // self-cumulate the subdivide result
            Subdivide_Combine(half_pos,half, half_pos,half, cache,
                sequencing.op_normal, sequencing.op_normal_flip,
                synth);
        }
        else
        {
            long part1 = half;
            long part2 = otherhalf>0?otherhalf:-otherhalf;

            size_t part1_pos = AssembleSequence_Subdivide(part1, cache, sequencing, synth);
            size_t part2_pos = AssembleSequence_Subdivide(part2, cache, sequencing, synth);

            FPO(fprintf(stderr, "Subdivide(%ld: %ld, %ld)\n", value, half, otherhalf));

            Subdivide_Combine(part1_pos,part1, part2_pos,part2, cache,
                otherhalf>0 ? sequencing.op_normal      : sequencing.op_inverse,
                otherhalf>0 ? sequencing.op_normal_flip : sequencing.op_inverse_flip,
                synth);
        }
        size_t stackpos = synth.GetStackTop()-1;
        cache.Remember(value, stackpos);
        cache.DumpContents();
        return stackpos;
    }

    void Subdivide_Combine(
        size_t apos, long aval,
        size_t bpos, long bval,
        PowiCache& cache,
        unsigned cumulation_opcode,
        unsigned cumulation_opcode_flip,
        ByteCodeSynth& synth)
    {
        /*FPO(fprintf(stderr, "== making result for (sp=%u, val=%d, needs=%d) and (sp=%u, val=%d, needs=%d), stacktop=%u\n",
            (unsigned)apos, aval, aval>=0 ? cache_needed[aval] : -1,
            (unsigned)bpos, bval, bval>=0 ? cache_needed[bval] : -1,
            (unsigned)synth.GetStackTop()));*/

        // Figure out whether we can trample a and b
        int a_needed = cache.UseGetNeeded(aval);
        int b_needed = cache.UseGetNeeded(bval);

        bool flipped = false;

        #define DUP_BOTH() do { \
            if(apos < bpos) { size_t tmp=apos; apos=bpos; bpos=tmp; flipped=!flipped; } \
            FPO(fprintf(stderr, "-> dup(%u) dup(%u) op\n", (unsigned)apos, (unsigned)bpos)); \
            synth.DoDup(apos); \
            synth.DoDup(apos==bpos ? synth.GetStackTop()-1 : bpos); } while(0)
        #define DUP_ONE(p) do { \
            FPO(fprintf(stderr, "-> dup(%u) op\n", (unsigned)p)); \
            synth.DoDup(p); \
        } while(0)

        if(a_needed > 0)
        {
            if(b_needed > 0)
            {
                // If they must both be preserved, make duplicates
                // First push the one that is at the larger stack
                // address. This increases the odds of possibly using cDup.
                DUP_BOTH();

                //SCENARIO 1:
                // Input:  x B A x x
                // Temp:   x B A x x A B
                // Output: x B A x x R
                //SCENARIO 2:
                // Input:  x A B x x
                // Temp:   x A B x x B A
                // Output: x A B x x R
            }
            else
            {
                // A must be preserved, but B can be trampled over

                // SCENARIO 1:
                //  Input:  x B x x A
                //   Temp:  x B x x A A B   (dup both, later first)
                //  Output: x B x x A R
                // SCENARIO 2:
                //  Input:  x A x x B
                //   Temp:  x A x x B A
                //  Output: x A x x R       -- only commutative cases
                // SCENARIO 3:
                //  Input:  x x x B A
                //   Temp:  x x x B A A B   (dup both, later first)
                //  Output: x x x B A R
                // SCENARIO 4:
                //  Input:  x x x A B
                //   Temp:  x x x A B A     -- only commutative cases
                //  Output: x x x A R
                // SCENARIO 5:
                //  Input:  x A B x x
                //   Temp:  x A B x x A B   (dup both, later first)
                //  Output: x A B x x R

                // if B is not at the top, dup both.
                if(bpos != synth.GetStackTop()-1)
                    DUP_BOTH();    // dup both
                else
                {
                    DUP_ONE(apos); // just dup A
                    flipped=!flipped;
                }
            }
        }
        else if(b_needed > 0)
        {
            // B must be preserved, but A can be trampled over
            // This is a mirror image of the a_needed>0 case, so I'll cut the chase
            if(apos != synth.GetStackTop()-1)
                DUP_BOTH();
            else
                DUP_ONE(bpos);
        }
        else
        {
            // Both can be trampled over.
            // SCENARIO 1:
            //  Input:  x B x x A
            //   Temp:  x B x x A B
            //  Output: x B x x R
            // SCENARIO 2:
            //  Input:  x A x x B
            //   Temp:  x A x x B A
            //  Output: x A x x R       -- only commutative cases
            // SCENARIO 3:
            //  Input:  x x x B A
            //  Output: x x x R         -- only commutative cases
            // SCENARIO 4:
            //  Input:  x x x A B
            //  Output: x x x R
            // SCENARIO 5:
            //  Input:  x A B x x
            //   Temp:  x A B x x A B   (dup both, later first)
            //  Output: x A B x x R
            // SCENARIO 6:
            //  Input:  x x x C
            //   Temp:  x x x C C   (c is both A and B)
            //  Output: x x x R

            if(apos == bpos && apos == synth.GetStackTop()-1)
                DUP_ONE(apos); // scenario 6
            else if(apos == synth.GetStackTop()-1 && bpos == synth.GetStackTop()-2)
            {
                FPO(fprintf(stderr, "-> op\n")); // scenario 3
                flipped=!flipped;
            }
            else if(apos == synth.GetStackTop()-2 && bpos == synth.GetStackTop()-1)
                FPO(fprintf(stderr, "-> op\n")); // scenario 4
            else if(apos == synth.GetStackTop()-1)
                DUP_ONE(bpos); // scenario 1
            else if(bpos == synth.GetStackTop()-1)
            {
                DUP_ONE(apos); // scenario 2
                flipped=!flipped;
            }
            else
                DUP_BOTH(); // scenario 5
        }
        // Add them together.
        synth.AddOperation(flipped ? cumulation_opcode_flip : cumulation_opcode, 2);
    }
}

namespace FPoptimizer_ByteCode
{
    void AssembleSequence(
        long count,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth)
    {
        if(count == 0)
            synth.PushImmed(sequencing.basevalue);
        else
        {
            if(count < 0)
            {
                synth.AddOperation(sequencing.op_flip, 1);
                count = -count;
            }

            if(count > 1)
            {
                /* To prevent calculating the same factors over and over again,
                 * we use a cache. */
                PowiCache cache;
                PlanNtimesCache(count, cache, 1);

                size_t stacktop_desired = synth.GetStackTop();

                cache.Start( synth.GetStackTop()-1 );

                FPO(fprintf(stderr, "Calculating result for %ld...\n", count));
                size_t res_stackpos = AssembleSequence_Subdivide(
                    count, cache, sequencing,
                    synth);

                size_t n_excess = synth.GetStackTop() - stacktop_desired;
                if(n_excess > 0 || res_stackpos != stacktop_desired-1)
                {
                    // Remove the cache values
                    synth.DoPopNMov(stacktop_desired-1, res_stackpos);
                }
            }
        }
    }
}
#include <cmath>
#include <list>
#include <algorithm>

#include <cmath> /* for CalculateResultBoundaries() */

#include "io/fparser/fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

//#define DEBUG_SUBSTITUTIONS

#ifdef DEBUG_SUBSTITUTIONS
namespace FPoptimizer_Grammar
{
    void DumpTree(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o = std::cout);
}
#endif

namespace FPoptimizer_CodeTree
{
    CodeTree::CodeTree()
        : RefCount(0), Opcode(), Params(), Hash(), Depth(1), Parent(), OptimizedUsing(0)
    {
    }

    CodeTree::CodeTree(double i)
        : RefCount(0), Opcode(cImmed), Params(), Hash(), Depth(1), Parent(), OptimizedUsing(0)
    {
        Value = i;
        Recalculate_Hash_NoRecursion();
    }

    CodeTree::~CodeTree()
    {
    }

    void CodeTree::Rehash(
        bool child_triggered)
    {
        /* If we were triggered by a parent, recurse to children */
        if(!child_triggered)
        {
            for(size_t a=0; a<Params.size(); ++a)
                Params[a].param->Rehash(false);
        }

        Recalculate_Hash_NoRecursion();

        /* If we were triggered by a child, recurse to the parent */
        if(child_triggered && Parent)
        {
            //assert(Parent->RefCount > 0);
            Parent->Rehash(true);
        }
    }

    struct ParamComparer
    {
        bool operator() (const CodeTree::Param& a, const CodeTree::Param& b) const
        {
            if(a.param->Depth != b.param->Depth)
                return a.param->Depth > b.param->Depth;
            if(a.sign != b.sign) return a.sign < b.sign;
            return a.param->Hash < b.param->Hash;
        }
    };

    void CodeTree::Sort()
    {
        /* If the tree is commutative, order the parameters
         * in a set order in order to make equality tests
         * efficient in the optimizer
         */
        switch(Opcode)
        {
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd:
            case cOr:
            case cEqual:
            case cNEqual:
                std::sort(Params.begin(), Params.end(), ParamComparer());
                break;
            case cLess:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cGreater; }
                break;
            case cLessOrEq:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cGreaterOrEq; }
                break;
            case cGreater:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cLess; }
                break;
            case cGreaterOrEq:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cLessOrEq; }
                break;
        }
    }

    void CodeTree::Sort_Recursive()
    {
        Sort();
        for(size_t a=0; a<Params.size(); ++a)
            Params[a].param->Sort_Recursive();
        Recalculate_Hash_NoRecursion();
    }

    void CodeTree::Recalculate_Hash_NoRecursion()
    {
        fphash_t NewHash = { Opcode * FPHASH_CONST(0x3A83A83A83A83A0),
                             Opcode * FPHASH_CONST(0x1131462E270012B)};
        Depth = 1;
        switch(Opcode)
        {
            case cImmed:
            {
                if(Value != 0.0)
                {
                    crc32_t crc = crc32::calc( (const unsigned char*) &Value,
                                                sizeof(Value) );
                    NewHash.hash1 ^= crc | (fphash_value_t(crc) << FPHASH_CONST(32));
                    NewHash.hash2 += ((~fphash_value_t(crc)) * 3) ^ 1234567;
                }
                break; // no params
            }
            case cVar:
                NewHash.hash1 ^= (Var<<24) | (Var>>24);
                NewHash.hash2 += (fphash_value_t(Var)*5) ^ 2345678;
                break; // no params
            case cFCall: case cPCall:
            {
                crc32_t crc = crc32::calc( (const unsigned char*) &Funcno, sizeof(Funcno) );
                NewHash.hash1 ^= (crc<<24) | (crc>>24);
                NewHash.hash2 += ((~fphash_value_t(crc)) * 7) ^ 3456789;
                /* passthru */
            }
            default:
            {
                size_t MaxChildDepth = 0;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    if(Params[a].param->Depth > MaxChildDepth)
                        MaxChildDepth = Params[a].param->Depth;

                    NewHash.hash1 += (1+Params[a].sign)*FPHASH_CONST(0x2492492492492492);
                    NewHash.hash1 *= FPHASH_CONST(1099511628211);
                    //assert(&*Params[a].param != this);
                    NewHash.hash1 += Params[a].param->Hash.hash1;

                    NewHash.hash2 += (3+Params[a].sign)*FPHASH_CONST(0x9ABCD801357);
                    NewHash.hash2 *= FPHASH_CONST(0xECADB912345);
                    NewHash.hash2 += (~Params[a].param->Hash.hash1) ^ 4567890;
                }
                Depth += MaxChildDepth;
            }
        }
        if(Hash != NewHash)
        {
            Hash = NewHash;
            OptimizedUsing = 0;
        }
    }

    CodeTree* CodeTree::Clone()
    {
        CodeTree* result = new CodeTree;
        result->Opcode = Opcode;
        switch(Opcode)
        {
            case cImmed:
                result->Value  = Value;
                break;
            case cVar:
                result->Var = Var;
                break;
            case cFCall: case cPCall:
                result->Funcno = Funcno;
                break;
        }
        result->SetParams(Params);
        result->Hash   = Hash;
        result->Depth  = Depth;
        //assert(Parent->RefCount > 0);
        result->Parent = Parent;
        return result;
    }

    void CodeTree::AddParam(const Param& param)
    {
        Params.push_back(param);
        Params.back().param->Parent = this;
    }

    void CodeTree::SetParams(const std::vector<Param>& RefParams, bool do_clone)
    {
        Params = RefParams;
        /**
        *** Note: The only reason we need to CLONE the children here
        ***       is because they must have the correct Parent field.
        ***       The Parent is required because of backward-recursive
        ***       hash regeneration. Is there any way around this?
        */

        for(size_t a=0; a<Params.size(); ++a)
        {
            if(do_clone) Params[a].param = Params[a].param->Clone();
            Params[a].param->Parent = this;
        }
    }

    void CodeTree::DelParam(size_t index)
    {
        Params.erase(Params.begin() + index);
    }

    CodeTree::MinMaxTree CodeTree::CalculateResultBoundaries() const
#ifdef DEBUG_SUBSTITUTIONS
    {
        MinMaxTree tmp = CalculateResultBoundaries_do();
        std::cout << std::flush;
        fprintf(stderr, "Estimated boundaries: %g%s .. %g%s: ",
            tmp.min, tmp.has_min?"":"(unknown)",
            tmp.max, tmp.has_max?"":"(unknown)");
        fflush(stderr);
        FPoptimizer_Grammar::DumpTree(*this);
        std::cout << std::flush;
        fprintf(stderr, " \n");
        fflush(stderr);
        return tmp;
    }
    CodeTree::MinMaxTree CodeTree::CalculateResultBoundaries_do() const
#endif
    {
        using namespace std;
        switch( (OPCODE) Opcode)
        {
            case cImmed:
                return MinMaxTree(Value, Value); // a definite value.
            case cEval:
            case cAnd:
            case cOr:
            case cNot:
            case cNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
            {
                /* These operations always produce truth values (0 or 1) */
                /* Narrowing them down is a matter of performing Constant optimization */
                return MinMaxTree( 0.0, 1.0 );
            }
            case cAbs:
            {
                /* cAbs always produces a positive value */
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min && m.has_max)
                {
                    if(m.min < 0.0 && m.max >= 0.0) // ex. -10..+6 or -6..+10
                    {
                        /* -x..+y: spans across zero. min=0, max=greater of |x| and |y|. */
                        double tmp = -m.min; if(tmp > m.max) m.max = tmp;
                        m.min = 0.0; m.has_min = true;
                    }
                    else if(m.min < 0.0) // ex. -10..-4
                        { double tmp = m.max; m.max = -m.min; m.min = -tmp; }
                }
                else if(!m.has_min && m.has_max && m.max < 0.0) // ex. -inf..-10
                {
                    m.min = fabs(m.max); m.has_min = true; m.has_max = false;
                }
                else if(!m.has_max && m.has_min && m.min > 0.0) // ex. +10..+inf
                {
                    m.min = fabs(m.min); m.has_min = true; m.has_max = false;
                }
                else // ex. -inf..+inf, -inf..+10, -10..+inf
                {
                    // all of these cover -inf..0, 0..+inf, or both
                    m.min = 0.0; m.has_min = true; m.has_max = false;
                }
                return m;
            }

            case cLog: /* Defined for 0.0 < x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) { if(m.min < 0.0) m.has_min = false; else m.min = log(m.min); } // No boundaries
                if(m.has_max) { if(m.max < 0.0) m.has_max = false; else m.max = log(m.max); }
                return m;
            }

            case cLog2: /* Defined for 0.0 < x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) { if(m.min < 0.0) m.has_min = false; else m.min = log(m.min)*CONSTANT_L2I; } // No boundaries
                if(m.has_max) { if(m.max < 0.0) m.has_max = false; else m.max = log(m.max)*CONSTANT_L2I; }
                return m;
            }

            case cAcosh: /* defined for             1.0 <  x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) { if(m.min <= 1.0) m.has_min = false; else m.min = fp_acosh(m.min); } // No boundaries
                if(m.has_max) { if(m.max <= 1.0) m.has_max = false; else m.max = fp_acosh(m.max); }
                return m;
            }
            case cAsinh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) m.min = fp_asinh(m.min); // No boundaries
                if(m.has_max) m.max = fp_asinh(m.max);
                return m;
            }
            case cAtanh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) m.min = fp_atanh(m.min); // No boundaries
                if(m.has_max) m.max = fp_atanh(m.max);
                return m;
            }
            case cAcos: /* defined for -1.0 <= x < 1, results within CONSTANT_PI..0 */
            {
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree( 0.0, CONSTANT_PI );
            }
            case cAsin: /* defined for -1.0 <= x < 1, results within -CONSTANT_PIHALF..CONSTANT_PIHALF */
            {
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree( -CONSTANT_PIHALF, CONSTANT_PIHALF );
            }
            case cAtan: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) m.min = atan(m.min); else { m.min = -CONSTANT_PIHALF; m.has_min = true; }
                if(m.has_max) m.max = atan(m.max); else { m.max =  CONSTANT_PIHALF; m.has_max = true; }
                return m;
            }
            case cAtan2: /* too complicated to estimate */
            {
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree(-CONSTANT_PI, CONSTANT_PI);
            }

            case cSin:
            case cCos:
            {
                /* Could be narrowed down from here,
                 * but it's too complicated due to
                 * the cyclic nature of the function. */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree(-1.0, 1.0);
            }
            case cTan:
            {
                /* Could be narrowed down from here,
                 * but it's too complicated due to
                 * the cyclic nature of the function */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree(); // (CONSTANT_NEG_INF, CONSTANT_POS_INF);
            }

            case cCeil:
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                m.max = std::ceil(m.max); // ceil() may increase the value, may not decrease
                return m;
            }
            case cFloor:
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                m.min = std::floor(m.min); // floor() may decrease the value, may not increase
                return m;
            }
            case cInt:
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                m.min = std::floor(m.min); // int() may either increase or decrease the value
                m.max = std::ceil(m.max); // for safety, we assume both
                return m;
            }
            case cSinh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) m.min = sinh(m.min); // No boundaries
                if(m.has_max) m.max = sinh(m.max);
                return m;
            }
            case cTanh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min) m.min = tanh(m.min); // No boundaries
                if(m.has_max) m.max = tanh(m.max);
                return m;
            }
            case cCosh: /* defined for all values -inf <= x <= inf, results within 1..inf */
            {
                MinMaxTree m = Params[0].param->CalculateResultBoundaries();
                if(m.has_min)
                {
                    if(m.has_max) // max, min
                    {
                        if(m.min >= 0.0 && m.max >= 0.0) // +x .. +y
                            { m.min = cosh(m.min); m.max = cosh(m.max); }
                        else if(m.min < 0.0 && m.max >= 0.0) // -x .. +y
                            { double tmp = cosh(m.min); m.max = cosh(m.max);
                              if(tmp > m.max) m.max = tmp;
                              m.min = 1.0; }
                        else // -x .. -y
                            { m.min = cosh(m.min); m.max = cosh(m.max);
                              std::swap(m.min, m.max); }
                    }
                    else // min, no max
                    {
                        if(m.min >= 0.0) // 0..inf -> 1..inf
                            { m.has_max = true; m.max = cosh(m.min); m.min = 1.0; }
                        else
                            { m.has_max = false; m.min = 1.0; } // Anything between 1..inf
                    }
                }
                else // no min
                {
                    m.has_min = true; m.min = 1.0; // always a lower boundary
                    if(m.has_max) // max, no min
                    {
                        m.min = cosh(m.max); // n..inf
                        m.has_max = false; // No upper boundary
                    }
                    else // no max, no min
                        m.has_max = false; // No upper boundary
                }
                return m;
            }

            case cIf:
            {
                // No guess which branch is chosen. Produce a spanning min & max.
                MinMaxTree res1 = Params[1].param->CalculateResultBoundaries();
                MinMaxTree res2 = Params[2].param->CalculateResultBoundaries();
                if(!res2.has_min) res1.has_min = false; else if(res2.min < res1.min) res1.min = res2.min;
                if(!res2.has_max) res1.has_max = false; else if(res2.max > res1.max) res1.max = res2.max;
                return res1;
            }

            case cMin:
            {
                bool has_unknown_min = false;
                bool has_unknown_max = false;

                MinMaxTree result;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    MinMaxTree m = Params[a].param->CalculateResultBoundaries();
                    if(!m.has_min)
                        has_unknown_min = true;
                    else if(!result.has_min || m.min < result.min)
                        result.min = m.min;

                    if(!m.has_max)
                        has_unknown_max = true;
                    else if(!result.has_max || m.max < result.max)
                        result.max = m.max;
                }
                if(has_unknown_min) result.has_min = false;
                if(has_unknown_max) result.has_max = false;
                return result;
            }
            case cMax:
            {
                bool has_unknown_min = false;
                bool has_unknown_max = false;

                MinMaxTree result;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    MinMaxTree m = Params[a].param->CalculateResultBoundaries();
                    if(!m.has_min)
                        has_unknown_min = true;
                    else if(!result.has_min || m.min > result.min)
                        result.min = m.min;

                    if(!m.has_max)
                        has_unknown_max = true;
                    else if(!result.has_max || m.max > result.max)
                        result.max = m.max;
                }
                if(has_unknown_min) result.has_min = false;
                if(has_unknown_max) result.has_max = false;
                return result;
            }
            case cAdd:
            {
                /* It's complicated. Follow the logic below. */
                /* Note: This also deals with the following opcodes:
                 *       cNeg, cSub, cRSub
                 */
                MinMaxTree result(0.0, 0.0);
                for(size_t a=0; a<Params.size(); ++a)
                {
                    const Param& p = Params[a];
                    MinMaxTree item = p.param->CalculateResultBoundaries();

                    if(Params[a].sign) return MinMaxTree(); // minus = unpredictable

                    if(item.has_min) result.min += item.min;
                    else             result.has_min = false;
                    if(item.has_max) result.max += item.max;
                    else             result.has_max = false;

                    if(!result.has_min && !result.has_max) break; // hopeless
                }
                if(result.has_min && result.has_max
                && result.min > result.max) std::swap(result.min, result.max);
                return result;
            }
            case cMul:
            {
                /* It's complicated. Follow the logic below. */
                /* Note: This also deals with the following opcodes:
                 *       cInv, cDiv, cRDiv, cRad, cDeg, cSqr
                 *       cCot, Sec, cCsc, cLog2, cLog10
                 */

                struct Value
                {
                    enum ValueType { Finite, MinusInf, PlusInf };
                    ValueType valueType;
                    double value;

                    Value(ValueType t): valueType(t), value(0) {}
                    Value(double v): valueType(Finite), value(v) {}

                    bool isNegative() const
                    {
                        return valueType == MinusInf ||
                            (valueType == Finite && value < 0.0);
                    }

                    void operator*=(const Value& rhs)
                    {
                        if(valueType == Finite && rhs.valueType == Finite)
                            value *= rhs.value;
                        else
                            valueType = (isNegative() != rhs.isNegative() ?
                                         MinusInf : PlusInf);
                    }

                    bool operator<(const Value& rhs) const
                    {
                        return
                            (valueType == MinusInf && rhs.valueType != MinusInf) ||
                            (valueType == Finite &&
                             (rhs.valueType == PlusInf ||
                              (rhs.valueType == Finite && value < rhs.value)));
                    }
                };

                struct MultiplicationRange
                {
                    Value minValue, maxValue;

                    MultiplicationRange():
                        minValue(Value::PlusInf),
                        maxValue(Value::MinusInf) {}

                    void multiply(Value value1, const Value& value2)
                    {
                        value1 *= value2;
                        if(value1 < minValue) minValue = value1;
                        if(maxValue < value1) maxValue = value1;
                    }
                };

                MinMaxTree result(1.0, 1.0);
                for(size_t a=0; a<Params.size(); ++a)
                {
                    const Param& p = Params[a];
                    MinMaxTree item = p.param->CalculateResultBoundaries();
                    if(!item.has_min && !item.has_max) return MinMaxTree(); // hopeless

                    if(Params[a].sign) return MinMaxTree(); // inversion = unpredictable

                    Value minValue0 = result.has_min ? Value(result.min) : Value(Value::MinusInf);
                    Value maxValue0 = result.has_max ? Value(result.max) : Value(Value::PlusInf);
                    Value minValue1 = item.has_min ? Value(item.min) : Value(Value::MinusInf);
                    Value maxValue1 = item.has_max ? Value(item.max) : Value(Value::PlusInf);

                    MultiplicationRange range;
                    range.multiply(minValue0, minValue1);
                    range.multiply(minValue0, maxValue1);
                    range.multiply(maxValue0, minValue1);
                    range.multiply(maxValue0, maxValue1);

                    if(range.minValue.valueType == Value::Finite)
                        result.min = range.minValue.value;
                    else result.has_min = false;

                    if(range.maxValue.valueType == Value::Finite)
                        result.max = range.maxValue.value;
                    else result.has_max = false;

                    if(!result.has_min && !result.has_max) break; // hopeless
                }
                if(result.has_min && result.has_max
                && result.min > result.max) std::swap(result.min, result.max);
                return result;
            }
            case cMod:
            {
                /* TODO: The boundaries of modulo operator could be estimated better. */

                MinMaxTree x = Params[0].param->CalculateResultBoundaries();
                MinMaxTree y = Params[1].param->CalculateResultBoundaries();

                if(y.has_max)
                {
                    if(y.max >= 0.0)
                    {
                        if(!x.has_min || x.min < 0)
                            return MinMaxTree(-y.max, y.max);
                        else
                            return MinMaxTree(0.0, y.max);
                    }
                    else
                    {
                        if(!x.has_max || x.max >= 0)
                            return MinMaxTree(y.max, -y.max);
                        else
                            return MinMaxTree(y.max, NEGATIVE_MAXIMUM);
                    }
                }
                else
                    return MinMaxTree();
            }
            case cPow:
            {
                if(Params[1].param->IsImmed() && FloatEqual(Params[1].param->GetImmed(), 0.0))
                {
                    // Note: This makes 0^0 evaluate into 1.
                    return MinMaxTree(1.0, 1.0); // x^0 = 1
                }
                if(Params[0].param->IsImmed() && FloatEqual(Params[0].param->GetImmed(), 0.0))
                {
                    // Note: This makes 0^0 evaluate into 0.
                    return MinMaxTree(0.0, 0.0); // 0^x = 0
                }
                if(Params[0].param->IsImmed() && FloatEqual(Params[0].param->GetImmed(), 1.0))
                {
                    return MinMaxTree(1.0, 1.0); // 1^x = 1
                }

                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                TriTruthValue p0_positivity =
                    (p0.has_min && p0.has_max)
                        ? ( (p0.min >= 0.0 && p0.max >= 0.0) ? IsAlways
                          : (p0.min <  0.0 && p0.max <  0.0) ? IsNever
                          : Unknown)
                        : Unknown;
                TriTruthValue p1_evenness = Params[1].param->GetEvennessInfo();

                /* If param0 IsAlways, the return value is also IsAlways */
                /* If param1 is even, the return value is IsAlways */
                /* If param1 is odd, the return value is same as param0's */
                /* If param0 is negative and param1 is not integer,
                 * the return value is imaginary (assumed Unknown)
                 *
                 * Illustrated in this truth table:
                 *  P=positive, N=negative
                 *  E=even, O=odd, U=not integer
                 *  *=unknown, X=invalid (unknown), x=maybe invalid (unknown)
                 *
                 *   param1: PE PO P* NE NO N* PU NU *
                 * param0:
                 *   PE      P  P  P  P  P  P  P  P  P
                 *   PO      P  P  P  P  P  P  P  P  P
                 *   PU      P  P  P  P  P  P  P  P  P
                 *   P*      P  P  P  P  P  P  P  P  P
                 *   NE      P  N  *  P  N  *  X  X  x
                 *   NO      P  N  *  P  N  *  X  X  x
                 *   NU      P  N  *  P  N  *  X  X  x
                 *   N*      P  N  *  P  N  *  X  X  x
                 *   *       P  *  *  P  *  *  x  x  *
                 *
                 * Note: This also deals with the following opcodes:
                 *       cSqrt  (param0, PU) (x^0.5)
                 *       cRSqrt (param0, NU) (x^-0.5)
                 *       cExp   (PU, param1) (CONSTANT_E^x)
                 */
                TriTruthValue result_positivity = Unknown;
                switch(p0_positivity)
                {
                    case IsAlways:
                        // e.g.   5^x = positive.
                        result_positivity = IsAlways;
                        break;
                    case IsNever:
                    {
                        result_positivity = p1_evenness;
                        break;
                    }
                    default:
                        switch(p1_evenness)
                        {
                            case IsAlways:
                                // e.g. x^( 4) = positive
                                // e.g. x^(-4) = positive
                                result_positivity = IsAlways;
                                break;
                            case IsNever:
                                break;
                            case Unknown:
                            {
                                /* If p1 is const non-integer,
                                 * assume the result is positive
                                 * though it may be NaN instead.
                                 */
                                if(Params[1].param->IsImmed()
                                && !Params[1].param->IsAlwaysInteger()
                                && Params[1].param->GetImmed() >= 0.0)
                                {
                                    result_positivity = IsAlways;
                                }
                                break;
                            }
                        }
                }
                switch(result_positivity)
                {
                    case IsAlways:
                    {
                        /* The result is always positive.
                         * Figure out whether we know the minimum value. */
                        double min = 0.0;
                        if(p0.has_min && p1.has_min)
                        {
                            min = pow(p0.min, p1.min);
                            if(p0.min < 0.0 && (!p1.has_max || p1.max >= 0.0) && min >= 0.0)
                                min = 0.0;
                        }
                        if(p0.has_min && p0.min >= 0.0 && p0.has_max && p1.has_max)
                        {
                            double max = pow(p0.max, p1.max);
                            if(min > max) std::swap(min, max);
                            return MinMaxTree(min, max);
                        }
                        return MinMaxTree(min, false);
                    }
                    case IsNever:
                    {
                        /* The result is always negative.
                         * TODO: Figure out whether we know the maximum value.
                         */
                        return MinMaxTree(false, NEGATIVE_MAXIMUM);
                    }
                    default:
                    {
                        /* It can be negative or positive.
                         * We know nothing about the boundaries. */
                        break;
                    }
                }
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.cc and thus
             * they will never occur in the calling context:
             */
            case cNeg: // converted into cAdd ~x
            case cInv: // converted into cMul ~x
            case cDiv: // converted into cMul ~x
            case cRDiv: // similar to above
            case cSub: // converted into cAdd ~x
            case cRSub: // similar to above
            case cRad: // converted into cMul x CONSTANT_RD
            case cDeg: // converted into cMul x CONSTANT_DR
            case cSqr: // converted into cMul x x
            case cExp: // converted into cPow CONSTANT_E x
            case cExp2: // converted into cPow 2 x
            case cSqrt: // converted into cPow x 0.5
            case cRSqrt: // converted into cPow x -0.5
            case cCot: // converted into cMul ~(cTan x)
            case cSec: // converted into cMul ~(cCos x)
            case cCsc: // converted into cMul ~(cSin x)
            case cLog10: // converted into cMul CONSTANT_L10I (cLog x)
                break; /* Should never occur */

            /* Opcodes that do not occur in the tree for other reasons */
            case cDup:
            case cFetch:
            case cPopNMov:
            case cNop:
            case cJump:
            case VarBegin:
                break; /* Should never occur */

            /* Opcodes that are completely unpredictable */
            case cVar:
            case cPCall:
            case cFCall:
#         ifndef FP_DISABLE_EVAL
            case cEval:
#endif
                break; // Cannot deduce


            //default:
                break;
        }
        return MinMaxTree(); /* Cannot deduce */
    }

    /* Is the value of this tree definitely odd(true) or even(false)? */
    CodeTree::TriTruthValue CodeTree::GetEvennessInfo() const
    {
        if(!IsImmed()) return Unknown;
        if(!IsLongIntegerImmed()) return Unknown;
        return (GetLongIntegerImmed() & 1) ? IsNever : IsAlways;
    }

    bool CodeTree::IsAlwaysInteger() const
    {
        switch( (OPCODE) Opcode)
        {
            case cImmed:
                return IsLongIntegerImmed();
            case cAnd:
            case cOr:
            case cNot:
            case cNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
                /* These operations always produce truth values (0 or 1) */
                return true; /* 0 and 1 are both integers */
            default:
                break;
        }
        return false; /* Don't know whether it's integer. */
    }

    bool CodeTree::IsAlwaysSigned(bool positive) const
    {
        MinMaxTree tmp = CalculateResultBoundaries();

        if(positive)
            return tmp.has_min && tmp.min >= 0.0
              && (!tmp.has_max || tmp.max >= 0.0);
        else
            return tmp.has_max && tmp.max < 0.0
              && (!tmp.has_min || tmp.min < 0.0);
    }

    bool CodeTree::IsIdenticalTo(const CodeTree& b) const
    {
        if(Hash   != b.Hash) return false; // a quick catch-all
        if(Opcode != b.Opcode) return false;
        switch(Opcode)
        {
            case cImmed: if(Value != b.Value) return false; return true;
            case cVar:   if(Var   != b.Var)   return false; return true;
            case cFCall:
            case cPCall: if(Funcno != b.Funcno) return false; break;
        }
        if(Params.size() != b.Params.size()) return false;
        for(size_t a=0; a<Params.size(); ++a)
        {
            if(Params[a].sign != b.Params[a].sign) return false;
            if(!Params[a].param->IsIdenticalTo(
               *b.Params[a].param)) return false;
        }
        return true;
    }

    bool    CodeTree::IsImmed() const { return Opcode == cImmed; }
    bool    CodeTree::IsVar()   const { return Opcode == cVar; }
}

#endif
/* This file is automatically generated. Do not edit... */
#include "io/fparser/fpconfig.h"
#include "io/fparser/fptypes.h"

using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

namespace
{
    const double clist[] =
    {
        2, /* 0 */
        -1, /* 1 */
        0.5, /* 2 */
        1, /* 3 */
        0.36787944117144233402427744294982403516769409179688, /* 4 */
        7.3890560989306504069418224389664828777313232421875, /* 5 */
        1.5707963267948965579989817342720925807952880859375, /* 6 */
        0, /* 7 */
        2.7182818284590450907955982984276488423347473144531, /* 8 */
        -2, /* 9 */
        0.13533528323661270231781372785917483270168304443359, /* 10 */
        0.017453292519943295474371680597869271878153085708618, /* 11 */
        57.29577951308232286464772187173366546630859375, /* 12 */
    };

    const ParamSpec plist[] =
    {
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 0 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 1 "x"	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 2    	*/
        {SubFunction , false, None  , 1, false, 0,	0 }, /* 3    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 4    	*/
        {SubFunction , false, None  , 1, false, 0,	1 }, /* 5    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 6    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 7 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	2 }, /* 8    	*/
        {SubFunction , false, None  , 1, false, 0,	3 }, /* 9    	*/
        {SubFunction , false, None  , 1, false, 0,	4 }, /* 10    	*/
        {SubFunction , false, None  , 1, false, 0,	0 }, /* 11    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 12    	*/
        {SubFunction , false, None  , 1, false, 0,	5 }, /* 13    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 14    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 15 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	6 }, /* 16    	*/
        {SubFunction , false, None  , 1, false, 0,	7 }, /* 17    	*/
        {SubFunction , false, None  , 1, false, 0,	8 }, /* 18    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 19    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 20 "x"	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 21    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 22 "x"	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 23    	*/
        {SubFunction , false, None  , 1, false, 0,	10 }, /* 24    	*/
        {SubFunction , false, None  , 1, false, 0,	11 }, /* 25    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 26    	*/
        {SubFunction , false, None  , 1, false, 0,	9 }, /* 27    	*/
        {SubFunction , false, None  , 1, false, 0,	12 }, /* 28    	*/
        {SubFunction , false, None  , 1, false, 0,	13 }, /* 29    	*/
        {SubFunction , false, None  , 1, false, 0,	14 }, /* 30    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 31    	*/
        {SubFunction , false, None  , 1, false, 0,	15 }, /* 32    	*/
        {NumConstant , false, None  , 1, false, 0,	4 }, /* 33    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 34 "x"	*/
        {NumConstant , false, None  , 1, false, 0,	5 }, /* 35    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 36 "x"	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 37    	*/
        {SubFunction , false, None  , 1, false, 0,	17 }, /* 38    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 39    	*/
        {SubFunction , false, None  , 1, false, 0,	16 }, /* 40    	*/
        {SubFunction , false, None  , 1, false, 0,	18 }, /* 41    	*/
        {SubFunction , false, None  , 1, false, 0,	19 }, /* 42    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 43    	*/
        {SubFunction , false, None  , 1, false, 0,	17 }, /* 44    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 45    	*/
        {SubFunction , false, None  , 1, false, 0,	16 }, /* 46    	*/
        {SubFunction , false, None  , 1, false, 0,	20 }, /* 47    	*/
        {SubFunction , false, None  , 1, false, 0,	21 }, /* 48    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 49    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 50    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 51    	*/
        {SubFunction , false, None  , 1, false, 0,	24 }, /* 52    	*/
        {SubFunction , false, None  , 1, false, 0,	25 }, /* 53    	*/
        {SubFunction , false, None  , 1, false, 0,	17 }, /* 54    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 55    	*/
        {SubFunction , false, None  , 1, false, 0,	17 }, /* 56    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 57    	*/
        {SubFunction , false, None  , 1, false, 0,	27 }, /* 58    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 59    	*/
        {SubFunction , false, None  , 1, false, 0,	26 }, /* 60    	*/
        {SubFunction , false, None  , 1, false, 0,	28 }, /* 61    	*/
        {SubFunction , false, None  , 1, false, 0,	29 }, /* 62    	*/
        {SubFunction , false, None  , 1, false, 0,	30 }, /* 63    	*/
        {SubFunction , false, None  , 1, false, 0,	31 }, /* 64    	*/
        {NamedHolder , false, None  , 1, false, Sign_Negative,	0 }, /* 65 "x"	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 66    	*/
        {SubFunction , false, None  , 1, false, 0,	32 }, /* 67    	*/
        {SubFunction , false, None  , 1, false, 0,	33 }, /* 68    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 69    	*/
        {SubFunction , false, None  , 1, false, 0,	34 }, /* 70    	*/
        {SubFunction , false, None  , 1, false, 0,	35 }, /* 71    	*/
        {NumConstant , false, Negate, 1, false, 0,	6 }, /* 72    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 73    	*/
        {SubFunction , false, None  , 1, false, 0,	36 }, /* 74    	*/
        {SubFunction , false, None  , 1, false, 0,	37 }, /* 75    	*/
        {SubFunction , false, None  , 1, false, 0,	38 }, /* 76    	*/
        {ImmedHolder , false, None  , 1, false, Sign_Negative,	0 }, /* 77    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 78    	*/
        {NumConstant , false, None  , 1, false, 0,	6 }, /* 79    	*/
        {SubFunction , false, None  , 1, false, 0,	39 }, /* 80    	*/
        {SubFunction , false, None  , 1, false, 0,	40 }, /* 81    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 82    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 83    	*/
        {SubFunction , false, None  , 1, false, 0,	41 }, /* 84    	*/
        {SubFunction , false, None  , 1, false, 0,	42 }, /* 85    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 86 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 87 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 88 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 89 "x"	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 90    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 91 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 92 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	43 }, /* 93    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 94 "y"	*/
        {NumConstant , false, None  , 1, false, 0,	7 }, /* 95    	*/
        {SubFunction , false, None  , 1, false, 0,	37 }, /* 96    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 97 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	44 }, /* 98    	*/
        {SubFunction , false, None  , 1, false, 0,	45 }, /* 99    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 100 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 101 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	46 }, /* 102    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 103 "y"	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 104    	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 105    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 106 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	48 }, /* 107    	*/
        {SubFunction , false, None  , 1, false, 0,	49 }, /* 108    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 109 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 110 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	50 }, /* 111    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 112 "y"	*/
        {NumConstant , false, None  , 1, false, 0,	7 }, /* 113    	*/
        {SubFunction , false, None  , 1, false, 0,	51 }, /* 114    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 115 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	52 }, /* 116    	*/
        {SubFunction , false, None  , 1, false, 0,	53 }, /* 117    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 118 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 119 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	54 }, /* 120    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 121 "y"	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 122    	*/
        {SubFunction , false, None  , 1, false, 0,	55 }, /* 123    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 124 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	56 }, /* 125    	*/
        {SubFunction , false, None  , 1, false, 0,	57 }, /* 126    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 127 "x"	*/
        {RestHolder  , false, None  , 1, false, 0,	3 }, /* 128    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 129 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	43 }, /* 130    	*/
        {SubFunction , false, None  , 1, false, 0,	58 }, /* 131    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 132 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	37 }, /* 133    	*/
        {SubFunction , false, None  , 1, false, 0,	59 }, /* 134    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 135 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	60 }, /* 136    	*/
        {SubFunction , false, None  , 1, false, 0,	61 }, /* 137    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 138 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	46 }, /* 139    	*/
        {SubFunction , false, None  , 1, false, 0,	62 }, /* 140    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 141 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 142    	*/
        {SubFunction , false, None  , 1, false, 0,	63 }, /* 143    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 144 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	64 }, /* 145    	*/
        {SubFunction , false, None  , 1, false, 0,	65 }, /* 146    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 147 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	50 }, /* 148    	*/
        {SubFunction , false, None  , 1, false, 0,	66 }, /* 149    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 150 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	51 }, /* 151    	*/
        {SubFunction , false, None  , 1, false, 0,	67 }, /* 152    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 153 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	68 }, /* 154    	*/
        {SubFunction , false, None  , 1, false, 0,	69 }, /* 155    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 156 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	54 }, /* 157    	*/
        {SubFunction , false, None  , 1, false, 0,	70 }, /* 158    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 159 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	55 }, /* 160    	*/
        {SubFunction , false, None  , 1, false, 0,	71 }, /* 161    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 162 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	72 }, /* 163    	*/
        {SubFunction , false, None  , 1, false, 0,	73 }, /* 164    	*/
        {SubFunction , false, None  , 1, false, 0,	74 }, /* 165    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 166 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 167 "z"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 168 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 169 "z"	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 170 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	75 }, /* 171    	*/
        {SubFunction , false, None  , 1, false, 0,	76 }, /* 172    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 173 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 174 "z"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 175 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 176 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 177 "z"	*/
        {SubFunction , false, None  , 1, false, 0,	77 }, /* 178    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 179 "x"	*/
        {NamedHolder , false, None  , 1, false, Value_EvenInt,	1 }, /* 180 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	78 }, /* 181    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 182 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	79 }, /* 183    	*/
        {SubFunction , false, None  , 1, false, 0,	80 }, /* 184    	*/
        {NamedHolder , false, None  , 1, false, Sign_Positive,	0 }, /* 185 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 186 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	81 }, /* 187    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 188 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 189    	*/
        {SubFunction , false, None  , 1, false, 0,	83 }, /* 190    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 191    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 192 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	84 }, /* 193    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 194    	*/
        {SubFunction , false, None  , 1, false, 0,	85 }, /* 195    	*/
        {SubFunction , false, None  , 1, false, 0,	86 }, /* 196    	*/
        {cLog        , false, Invert, 1, false, 1,	191 }, /* 197    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 198 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	87 }, /* 199    	*/
        {cLog        , false, None  , 1, false, 1,	191 }, /* 200    	*/
        {SubFunction , false, None  , 1, false, 0,	88 }, /* 201    	*/
        {SubFunction , false, None  , 1, false, 0,	89 }, /* 202    	*/
        {ImmedHolder , false, None  , 1, false, Sign_Positive,	0 }, /* 203    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 204    	*/
        {SubFunction , false, None  , 1, false, 0,	90 }, /* 205    	*/
        {SubFunction , false, None  , 1, false, 0,	86 }, /* 206    	*/
        {cLog        , false, None  , 1, false, 1,	191 }, /* 207    	*/
        {SubFunction , false, None  , 1, false, 0,	91 }, /* 208    	*/
        {NumConstant , false, None  , 1, false, 0,	8 }, /* 209    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 210    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 211    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 212    	*/
        {NumConstant , false, None  , 1, false, 0,	8 }, /* 213    	*/
        {SubFunction , false, None  , 1, false, 0,	92 }, /* 214    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 215 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 216    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 217    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 218    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 219 "x"	*/
        {cLog        , false, None  , 1, false, 1,	191 }, /* 220    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 221    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 222    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 223    	*/
        {SubFunction , false, None  , 1, false, 0,	93 }, /* 224    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 225    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 226    	*/
        {cPow        , false, None  , 1, false, 2,	225 }, /* 227    	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 228    	*/
        {cLog        , false, Invert, 1, false, 1,	191 }, /* 229    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 230    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 231    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 232    	*/
        {SubFunction , false, None  , 1, false, 0,	94 }, /* 233    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 234    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 235    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 236    	*/
        {SubFunction , false, None  , 1, false, 0,	95 }, /* 237    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 238 "x"	*/
        {cExp        , false, None  , 1, false, 1,	191 }, /* 239    	*/
        {SubFunction , false, None  , 1, false, 0,	96 }, /* 240    	*/
        {cLog        , false, None  , 1, false, 1,	221 }, /* 241    	*/
        {SubFunction , false, None  , 1, false, 0,	97 }, /* 242    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 243    	*/
        {SubFunction , false, None  , 1, false, 0,	98 }, /* 244    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 245    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 246    	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 247 "z"	*/
        {SubFunction , false, None  , 1, false, 0,	99 }, /* 248    	*/
        {SubFunction , false, None  , 1, false, 0,	33 }, /* 249    	*/
        {NamedHolder , false, None  , 1, false, Value_EvenInt,	1 }, /* 250 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 251 "x"	*/
        {NamedHolder , false, None  , 1, false, Value_OddInt,	1 }, /* 252 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	100 }, /* 253    	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 254 "z"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 255 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	101 }, /* 256    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 257 "x"	*/
        {NamedHolder , false, None  , 1, false, Value_NonInteger,	1 }, /* 258 "y"	*/
        {SubFunction , false, None  , 1, false, 0,	102 }, /* 259    	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 260 "z"	*/
        {SubFunction , false, None  , 1, false, 0,	81 }, /* 261    	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 262 "z"	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 263    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 264    	*/
        {SubFunction , false, None  , 1, false, 0,	103 }, /* 265    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 266    	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 267    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 268    	*/
        {cPow        , false, None  , 1, false, 2,	225 }, /* 269    	*/
        {SubFunction , false, None  , 1, false, 0,	104 }, /* 270    	*/
        {SubFunction , false, None  , 1, false, 0,	105 }, /* 271    	*/
        {NamedHolder , false, None  , 1, false, Sign_NoIdea,	0 }, /* 272 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 273 "y"	*/
        {SubFunction , false, None  , 1, false, 8,	106 }, /* 274    	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 275 "z"	*/
        {SubFunction , false, None  , 1, false, 0,	33 }, /* 276    	*/
        {SubFunction , false, None  , 1, false, 0,	101 }, /* 277    	*/
        {SubFunction , false, None  , 1, false, 0,	107 }, /* 278    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 279    	*/
        {SubFunction , false, None  , 1, false, 0,	108 }, /* 280    	*/
        {SubFunction , false, None  , 1, false, 0,	109 }, /* 281    	*/
        {SubFunction , false, None  , 1, false, 0,	110 }, /* 282    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 283    	*/
        {SubFunction , false, None  , 1, false, 0,	42 }, /* 284    	*/
        {SubFunction , false, None  , 1, false, 0,	111 }, /* 285    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 286    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 287    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 288    	*/
        {SubFunction , false, None  , 1, false, 0,	112 }, /* 289    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 290    	*/
        {SubFunction , false, None  , 1, false, 0,	113 }, /* 291    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 292    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 293    	*/
        {SubFunction , false, None  , 1, false, 0,	114 }, /* 294    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 295    	*/
        {SubFunction , false, None  , 1, false, 0,	114 }, /* 296    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 297    	*/
        {SubFunction , false, None  , 1, false, 0,	115 }, /* 298    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 299    	*/
        {SubFunction , false, None  , 1, false, 0,	116 }, /* 300    	*/
        {SubFunction , false, None  , 1, false, 0,	117 }, /* 301    	*/
        {SubFunction , false, None  , 1, false, 0,	118 }, /* 302    	*/
        {SubFunction , false, None  , 1, false, 0,	114 }, /* 303    	*/
        {SubFunction , false, None  , 1, false, 0,	112 }, /* 304    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 305 "x"	*/
        {RestHolder  , false, None  , 1, false, 0,	5 }, /* 306    	*/
        {SubFunction , false, None  , 1, false, 0,	62 }, /* 307    	*/
        {SubFunction , false, None  , 1, false, 0,	119 }, /* 308    	*/
        {SubFunction , false, None  , 1, false, 0,	63 }, /* 309    	*/
        {SubFunction , false, None  , 1, false, 0,	120 }, /* 310    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 311 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	121 }, /* 312    	*/
        {SubFunction , false, None  , 1, false, 0,	122 }, /* 313    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 314 "x"	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 315    	*/
        {SubFunction , false, None  , 1, false, 0,	123 }, /* 316    	*/
        {RestHolder  , false, None  , 1, false, 0,	3 }, /* 317    	*/
        {SubFunction , false, None  , 1, false, 0,	124 }, /* 318    	*/
        {SubFunction , false, None  , 1, false, 0,	119 }, /* 319    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 320    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 321    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 322 "x"	*/
        {cMin        , false, None  , 1, false, 2,	320 }, /* 323    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 324    	*/
        {cMin        , false, Negate, 1, false, 2,	320 }, /* 325    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 326 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	126 }, /* 327    	*/
        {RestHolder  , false, None  , 1, false, 0,	3 }, /* 328    	*/
        {SubFunction , false, None  , 1, false, 0,	127 }, /* 329    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 330    	*/
        {cMin        , false, Negate, 1, false, 2,	320 }, /* 331    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 332 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	129 }, /* 333    	*/
        {RestHolder  , false, None  , 1, false, 0,	5 }, /* 334    	*/
        {SubFunction , false, None  , 1, false, 0,	130 }, /* 335    	*/
        {SubFunction , false, None  , 1, false, 0,	128 }, /* 336    	*/
        {SubFunction , false, None  , 1, false, 0,	131 }, /* 337    	*/
        {SubFunction , false, None  , 1, false, 0,	125 }, /* 338    	*/
        {SubFunction , false, None  , 1, false, 0,	132 }, /* 339    	*/
        {SubFunction , false, None  , 1, false, 0,	133 }, /* 340    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 341 "x"	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 342    	*/
        {SubFunction , false, None  , 1, false, 0,	134 }, /* 343    	*/
        {RestHolder  , false, None  , 1, false, 0,	5 }, /* 344    	*/
        {SubFunction , false, None  , 1, false, 0,	124 }, /* 345    	*/
        {SubFunction , false, None  , 1, false, 0,	135 }, /* 346    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 347 "x"	*/
        {cMin        , false, None  , 1, false, 2,	225 }, /* 348    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 349    	*/
        {cMin        , false, Negate, 1, false, 2,	225 }, /* 350    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 351 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	137 }, /* 352    	*/
        {RestHolder  , false, None  , 1, false, 0,	3 }, /* 353    	*/
        {SubFunction , false, None  , 1, false, 0,	138 }, /* 354    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 355    	*/
        {cMin        , false, Negate, 1, false, 2,	225 }, /* 356    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 357 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	140 }, /* 358    	*/
        {RestHolder  , false, None  , 1, false, 0,	5 }, /* 359    	*/
        {SubFunction , false, None  , 1, false, 0,	141 }, /* 360    	*/
        {SubFunction , false, None  , 1, false, 0,	139 }, /* 361    	*/
        {SubFunction , false, None  , 1, false, 0,	142 }, /* 362    	*/
        {SubFunction , false, None  , 1, false, 0,	136 }, /* 363    	*/
        {SubFunction , false, None  , 1, false, 0,	143 }, /* 364    	*/
        {SubFunction , false, None  , 1, false, 0,	144 }, /* 365    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 366    	*/
        {SubFunction , false, None  , 1, false, 0,	145 }, /* 367    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 368    	*/
        {SubFunction , false, None  , 1, false, 0,	147 }, /* 369    	*/
        {SubFunction , false, None  , 1, false, 0,	146 }, /* 370    	*/
        {SubFunction , false, None  , 1, false, 0,	148 }, /* 371    	*/
        {SubFunction , false, None  , 1, false, 0,	149 }, /* 372    	*/
        {SubFunction , false, None  , 1, false, 0,	150 }, /* 373    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 374    	*/
        {SubFunction , false, None  , 1, false, 0,	147 }, /* 375    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 376    	*/
        {SubFunction , false, None  , 1, false, 0,	146 }, /* 377    	*/
        {SubFunction , false, None  , 1, false, 0,	151 }, /* 378    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 379    	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 380 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 381 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	152 }, /* 382    	*/
        {SubFunction , false, None  , 1, false, 0,	153 }, /* 383    	*/
        {SubFunction , false, None  , 1, false, 0,	154 }, /* 384    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 385    	*/
        {SubFunction , false, None  , 1, false, 0,	145 }, /* 386    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 387    	*/
        {SubFunction , false, None  , 1, false, 0,	147 }, /* 388    	*/
        {SubFunction , false, None  , 1, false, 0,	155 }, /* 389    	*/
        {SubFunction , false, None  , 1, false, 0,	156 }, /* 390    	*/
        {SubFunction , false, None  , 1, false, 0,	157 }, /* 391    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 392    	*/
        {SubFunction , false, None  , 1, false, 0,	147 }, /* 393    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 394    	*/
        {SubFunction , false, None  , 1, false, 0,	155 }, /* 395    	*/
        {SubFunction , false, None  , 1, false, 0,	158 }, /* 396    	*/
        {SubFunction , false, None  , 1, false, 0,	159 }, /* 397    	*/
        {ImmedHolder , false, None  , 1, false, Oneness_NotOne,	0 }, /* 398    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 399 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	160 }, /* 400    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 401    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 402    	*/
        {cAbs        , false, Invert, 1, false, 1,	191 }, /* 403    	*/
        {cMul        , false, None  , 1, false, 2,	402 }, /* 404    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 405 "x"	*/
        {cMul        , false, Negate, 1, false, 2,	402 }, /* 406    	*/
        {SubFunction , false, None  , 1, false, 0,	161 }, /* 407    	*/
        {cAbs        , false, None  , 1, false, 1,	191 }, /* 408    	*/
        {SubFunction , false, None  , 1, false, 0,	162 }, /* 409    	*/
        {SubFunction , false, None  , 1, false, 0,	163 }, /* 410    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 411    	*/
        {SubFunction , false, None  , 1, false, 0,	145 }, /* 412    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 413    	*/
        {SubFunction , false, None  , 1, false, 0,	164 }, /* 414    	*/
        {SubFunction , false, None  , 1, false, 0,	148 }, /* 415    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 416    	*/
        {SubFunction , false, None  , 1, false, 0,	154 }, /* 417    	*/
        {SubFunction , false, None  , 1, false, 0,	165 }, /* 418    	*/
        {SubFunction , false, None  , 1, false, 0,	164 }, /* 419    	*/
        {SubFunction , false, None  , 1, false, 0,	151 }, /* 420    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 421    	*/
        {SubFunction , false, None  , 1, false, 0,	150 }, /* 422    	*/
        {SubFunction , false, None  , 1, false, 0,	166 }, /* 423    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 424    	*/
        {SubFunction , false, None  , 1, false, 0,	145 }, /* 425    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 426    	*/
        {SubFunction , false, None  , 1, false, 0,	167 }, /* 427    	*/
        {SubFunction , false, None  , 1, false, 0,	158 }, /* 428    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 429    	*/
        {SubFunction , false, None  , 1, false, 0,	157 }, /* 430    	*/
        {SubFunction , false, None  , 1, false, 0,	168 }, /* 431    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 432    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 433 "x"	*/
        {ImmedHolder , false, None  , 1, false, Oneness_NotOne,	0 }, /* 434    	*/
        {SubFunction , false, None  , 1, false, 0,	169 }, /* 435    	*/
        {NamedHolder , false, None  , 1, true , 0,	0 }, /* 436 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	46 }, /* 437    	*/
        {NamedHolder , false, None  , 1, true , 0,	0 }, /* 438 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 439    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 440 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	170 }, /* 441    	*/
        {SubFunction , false, None  , 1, false, 0,	171 }, /* 442    	*/
        {NamedHolder , false, None  , 2, true , 0,	0 }, /* 443 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 444 "x"	*/
        {NamedHolder , false, None  , 2, true , 0,	0 }, /* 445 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	172 }, /* 446    	*/
        {NamedHolder , false, None  , 1, false, 0,	3 }, /* 447 "a"	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 448    	*/
        {NamedHolder , false, None  , 1, false, 0,	4 }, /* 449 "b"	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 450    	*/
        {NamedHolder , false, None  , 1, false, 0,	3 }, /* 451 "a"	*/
        {NamedHolder , false, None  , 1, false, 0,	4 }, /* 452 "b"	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 453    	*/
        {SubFunction , false, None  , 1, false, 0,	173 }, /* 454    	*/
        {SubFunction , false, None  , 1, false, 0,	174 }, /* 455    	*/
        {SubFunction , false, None  , 1, false, 0,	175 }, /* 456    	*/
        {SubFunction , false, None  , 1, false, 0,	176 }, /* 457    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 458    	*/
        {SubFunction , false, None  , 1, false, 0,	47 }, /* 459    	*/
        {NumConstant , false, None  , 1, false, 0,	9 }, /* 460    	*/
        {NamedHolder , false, None  , 1, false, 0,	3 }, /* 461 "a"	*/
        {NamedHolder , false, None  , 1, false, 0,	4 }, /* 462 "b"	*/
        {SubFunction , false, None  , 1, false, 0,	178 }, /* 463    	*/
        {SubFunction , false, None  , 1, false, 0,	177 }, /* 464    	*/
        {SubFunction , false, None  , 1, false, 0,	179 }, /* 465    	*/
        {SubFunction , false, None  , 1, false, 0,	82 }, /* 466    	*/
        {ImmedHolder , false, Invert, 1, false, 0,	0 }, /* 467    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 468    	*/
        {SubFunction , false, None  , 1, false, 0,	180 }, /* 469    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 470    	*/
        {SubFunction , false, None  , 1, false, 0,	181 }, /* 471    	*/
        {SubFunction , false, None  , 1, false, 0,	182 }, /* 472    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 473 "x"	*/
        {RestHolder  , false, None  , 1, false, 0,	2 }, /* 474    	*/
        {SubFunction , false, None  , 1, false, 0,	183 }, /* 475    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 476    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 477 "x"	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 478    	*/
        {SubFunction , false, None  , 1, false, 0,	184 }, /* 479    	*/
        {SubFunction , false, None  , 1, false, 0,	185 }, /* 480    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 481    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 482    	*/
        {SubFunction , false, None  , 1, false, 0,	186 }, /* 483    	*/
        {SubFunction , false, None  , 1, false, 0,	187 }, /* 484    	*/
        {SubFunction , false, None  , 1, false, 0,	188 }, /* 485    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 486 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	189 }, /* 487    	*/
        {SubFunction , false, None  , 1, false, 0,	190 }, /* 488    	*/
        {SubFunction , false, None  , 1, false, 0,	188 }, /* 489    	*/
        {SubFunction , false, None  , 1, false, 0,	191 }, /* 490    	*/
        {SubFunction , false, None  , 1, false, 0,	117 }, /* 491    	*/
        {NamedHolder , false, None  , 1, false, 0,	2 }, /* 492 "z"	*/
        {SubFunction , false, None  , 1, false, 0,	192 }, /* 493    	*/
        {ImmedHolder , false, None  , 1, false, Oneness_NotOne,	0 }, /* 494    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 495    	*/
        {SubFunction , false, None  , 1, false, 0,	193 }, /* 496    	*/
        {RestHolder  , false, None  , 1, false, 0,	3 }, /* 497    	*/
        {SubFunction , false, None  , 1, false, 0,	194 }, /* 498    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 499    	*/
        {cMul        , false, None  , 1, false, 2,	225 }, /* 500    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 501    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 502    	*/
        {SubFunction , false, None  , 1, false, 0,	59 }, /* 503    	*/
        {SubFunction , false, None  , 1, false, 0,	195 }, /* 504    	*/
        {SubFunction , false, None  , 1, false, 0,	196 }, /* 505    	*/
        {SubFunction , false, None  , 1, false, 0,	197 }, /* 506    	*/
        {SubFunction , false, None  , 1, false, 0,	198 }, /* 507    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 508    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 509    	*/
        {SubFunction , false, None  , 1, false, 0,	37 }, /* 510    	*/
        {cMul        , false, None  , 1, false, 2,	225 }, /* 511    	*/
        {SubFunction , false, None  , 1, false, 0,	199 }, /* 512    	*/
        {SubFunction , false, None  , 1, false, 0,	200 }, /* 513    	*/
        {NamedHolder , false, None  , 1, true , 0,	0 }, /* 514 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	187 }, /* 515    	*/
        {NamedHolder , false, None  , 1, true , 0,	0 }, /* 516 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	1 }, /* 517 "y"	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 518 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	201 }, /* 519    	*/
        {SubFunction , false, None  , 1, false, 0,	202 }, /* 520    	*/
        {SubFunction , false, None  , 1, false, 0,	203 }, /* 521    	*/
        {SubFunction , false, None  , 1, false, 0,	204 }, /* 522    	*/
        {SubFunction , false, None  , 1, false, 0,	205 }, /* 523    	*/
        {SubFunction , false, None  , 1, false, 0,	206 }, /* 524    	*/
        {SubFunction , false, None  , 1, false, 0,	207 }, /* 525    	*/
        {SubFunction , false, None  , 1, false, 0,	208 }, /* 526    	*/
        {SubFunction , false, None  , 1, false, 0,	209 }, /* 527    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 528    	*/
        {RestHolder  , true , None  , 1, false, 0,	2 }, /* 529    	*/
        {SubFunction , false, None  , 1, false, 0,	210 }, /* 530    	*/
        {RestHolder  , true , None  , 1, false, 0,	1 }, /* 531    	*/
        {RestHolder  , false, None  , 1, false, 0,	2 }, /* 532    	*/
        {SubFunction , false, None  , 1, false, 0,	211 }, /* 533    	*/
        {SubFunction , false, None  , 1, false, 0,	212 }, /* 534    	*/
        {SubFunction , false, None  , 1, false, 0,	213 }, /* 535    	*/
        {RestHolder  , false, None  , 1, false, 0,	2 }, /* 536    	*/
        {RestHolder  , true , None  , 1, false, 0,	1 }, /* 537    	*/
        {SubFunction , true , None  , 1, false, 0,	214 }, /* 538    	*/
        {NamedHolder , true , None  , 1, false, 0,	0 }, /* 539 "x"	*/
        {SubFunction , true , None  , 1, false, 0,	204 }, /* 540    	*/
        {SubFunction , true , None  , 1, false, 0,	205 }, /* 541    	*/
        {SubFunction , true , None  , 1, false, 0,	206 }, /* 542    	*/
        {SubFunction , true , None  , 1, false, 0,	208 }, /* 543    	*/
        {SubFunction , true , None  , 1, false, 0,	209 }, /* 544    	*/
        {SubFunction , true , None  , 1, false, 0,	207 }, /* 545    	*/
        {SubFunction , true , None  , 1, false, 0,	30 }, /* 546    	*/
        {SubFunction , true , None  , 1, false, 0,	31 }, /* 547    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 548 "x"	*/
        {NamedHolder , true , None  , 1, false, 0,	0 }, /* 549 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	206 }, /* 550    	*/
        {SubFunction , false, None  , 1, false, 0,	209 }, /* 551    	*/
        {SubFunction , false, None  , 1, false, 0,	208 }, /* 552    	*/
        {SubFunction , false, None  , 1, false, 0,	206 }, /* 553    	*/
        {SubFunction , false, None  , 1, false, 0,	208 }, /* 554    	*/
        {SubFunction , false, None  , 1, false, 0,	207 }, /* 555    	*/
        {SubFunction , false, None  , 1, false, 0,	207 }, /* 556    	*/
        {SubFunction , false, None  , 1, false, 0,	209 }, /* 557    	*/
        {NamedHolder , true , None  , 1, false, 0,	0 }, /* 558 "x"	*/
        {NamedHolder , true , None  , 1, false, 0,	0 }, /* 559 "x"	*/
        {NamedHolder , false, None  , 1, false, 0,	4 }, /* 560 "b"	*/
        {NamedHolder , false, None  , 1, false, 0,	5 }, /* 561 "c"	*/
        {NamedHolder , false, None  , 1, false, 0,	3 }, /* 562 "a"	*/
        {NamedHolder , false, None  , 1, false, 0,	5 }, /* 563 "c"	*/
        {SubFunction , false, None  , 1, false, 0,	204 }, /* 564    	*/
        {SubFunction , false, None  , 1, false, 0,	215 }, /* 565    	*/
        {SubFunction , false, None  , 1, false, 0,	216 }, /* 566    	*/
        {SubFunction , true , None  , 1, false, 0,	217 }, /* 567    	*/
        {SubFunction , false, None  , 1, false, 0,	206 }, /* 568    	*/
        {SubFunction , false, None  , 1, false, 0,	204 }, /* 569    	*/
        {SubFunction , false, None  , 1, false, 0,	208 }, /* 570    	*/
        {SubFunction , false, None  , 1, false, 0,	204 }, /* 571    	*/
        {SubFunction , false, None  , 1, false, 0,	209 }, /* 572    	*/
        {SubFunction , false, None  , 1, false, 0,	204 }, /* 573    	*/
        {SubFunction , false, None  , 1, false, 0,	207 }, /* 574    	*/
        {SubFunction , false, None  , 1, false, 0,	204 }, /* 575    	*/
        {SubFunction , false, None  , 1, false, 0,	218 }, /* 576    	*/
        {SubFunction , false, None  , 1, false, 0,	219 }, /* 577    	*/
        {SubFunction , false, None  , 1, false, 0,	220 }, /* 578    	*/
        {SubFunction , false, None  , 1, false, 0,	221 }, /* 579    	*/
        {SubFunction , false, None  , 1, false, 0,	222 }, /* 580    	*/
        {SubFunction , false, None  , 1, false, 0,	223 }, /* 581    	*/
        {SubFunction , false, None  , 1, false, 0,	224 }, /* 582    	*/
        {SubFunction , false, None  , 1, false, 0,	225 }, /* 583    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 584    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 585 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	226 }, /* 586    	*/
        {SubFunction , false, None  , 1, false, 0,	227 }, /* 587    	*/
        {SubFunction , false, None  , 1, false, 0,	228 }, /* 588    	*/
        {SubFunction , false, None  , 1, false, 0,	229 }, /* 589    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 590    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 591 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	230 }, /* 592    	*/
        {SubFunction , false, None  , 1, false, 0,	231 }, /* 593    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 594    	*/
        {SubFunction , false, None  , 1, false, 0,	233 }, /* 595    	*/
        {SubFunction , false, None  , 1, false, 0,	234 }, /* 596    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 597    	*/
        {SubFunction , false, None  , 1, false, 0,	232 }, /* 598    	*/
        {SubFunction , false, None  , 1, false, 0,	235 }, /* 599    	*/
        {SubFunction , false, None  , 1, false, 0,	236 }, /* 600    	*/
        {SubFunction , false, None  , 1, false, 0,	237 }, /* 601    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 602    	*/
        {SubFunction , false, None  , 1, false, 0,	238 }, /* 603    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 604    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 605 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	239 }, /* 606    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 607    	*/
        {SubFunction , false, None  , 1, false, 0,	169 }, /* 608    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 609    	*/
        {SubFunction , false, None  , 1, false, 0,	241 }, /* 610    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 611    	*/
        {SubFunction , false, None  , 1, false, 0,	240 }, /* 612    	*/
        {SubFunction , false, None  , 1, false, 0,	242 }, /* 613    	*/
        {SubFunction , false, None  , 1, false, 0,	243 }, /* 614    	*/
        {SubFunction , false, None  , 1, false, 0,	244 }, /* 615    	*/
        {SubFunction , false, None  , 1, false, 0,	245 }, /* 616    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 617    	*/
        {SubFunction , false, None  , 1, false, 0,	246 }, /* 618    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 619    	*/
        {SubFunction , false, None  , 1, false, 0,	247 }, /* 620    	*/
        {SubFunction , false, None  , 1, false, 0,	248 }, /* 621    	*/
        {SubFunction , false, None  , 1, false, 0,	249 }, /* 622    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 623    	*/
        {SubFunction , false, None  , 1, false, 0,	250 }, /* 624    	*/
        {SubFunction , false, None  , 1, false, 0,	251 }, /* 625    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 626    	*/
        {SubFunction , false, None  , 1, false, 0,	252 }, /* 627    	*/
        {SubFunction , false, None  , 1, false, 0,	253 }, /* 628    	*/
        {NumConstant , false, None  , 1, false, 0,	8 }, /* 629    	*/
        {SubFunction , false, None  , 1, false, 0,	228 }, /* 630    	*/
        {SubFunction , false, None  , 1, false, 0,	254 }, /* 631    	*/
        {NumConstant , false, None  , 1, false, 0,	8 }, /* 632    	*/
        {SubFunction , false, None  , 1, false, 0,	223 }, /* 633    	*/
        {SubFunction , false, None  , 1, false, 0,	255 }, /* 634    	*/
        {SubFunction , false, None  , 1, false, 0,	256 }, /* 635    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 636    	*/
        {cLog        , false, None  , 1, false, 1,	191 }, /* 637    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 638    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 639 "x"	*/
        {cMul        , false, None  , 1, false, 2,	637 }, /* 640    	*/
        {SubFunction , false, None  , 1, false, 0,	257 }, /* 641    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 642    	*/
        {NumConstant , false, None  , 1, false, 0,	2 }, /* 643    	*/
        {cPow        , false, None  , 1, false, 2,	642 }, /* 644    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 645 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	258 }, /* 646    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 647    	*/
        {SubFunction , false, None  , 1, false, 0,	259 }, /* 648    	*/
        {SubFunction , false, None  , 1, false, 0,	260 }, /* 649    	*/
        {SubFunction , false, None  , 1, false, 0,	256 }, /* 650    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 651    	*/
        {SubFunction , false, None  , 1, false, 0,	261 }, /* 652    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 653    	*/
        {SubFunction , false, None  , 1, false, 0,	259 }, /* 654    	*/
        {SubFunction , false, None  , 1, false, 0,	262 }, /* 655    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 656    	*/
        {SubFunction , false, None  , 1, false, 0,	264 }, /* 657    	*/
        {NumConstant , false, None  , 1, false, 0,	8 }, /* 658    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 659 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	265 }, /* 660    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 661    	*/
        {SubFunction , false, None  , 1, false, 0,	265 }, /* 662    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 663    	*/
        {SubFunction , false, None  , 1, false, 0,	266 }, /* 664    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 665    	*/
        {SubFunction , false, None  , 1, false, 0,	264 }, /* 666    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 667    	*/
        {SubFunction , false, None  , 1, false, 0,	267 }, /* 668    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 669    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 670 "x"	*/
        {ImmedHolder , false, Invert, 1, false, 0,	1 }, /* 671    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 672 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	268 }, /* 673    	*/
        {SubFunction , false, None  , 1, false, 0,	269 }, /* 674    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 675 "x"	*/
        {cLog        , false, None  , 1, false, 1,	221 }, /* 676    	*/
        {SubFunction , false, None  , 1, false, 0,	270 }, /* 677    	*/
        {SubFunction , false, None  , 1, false, 0,	271 }, /* 678    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 679    	*/
        {SubFunction , false, None  , 1, false, 0,	272 }, /* 680    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 681    	*/
        {SubFunction , false, None  , 1, false, 0,	269 }, /* 682    	*/
        {SubFunction , false, None  , 1, false, 0,	268 }, /* 683    	*/
        {SubFunction , false, None  , 1, false, 0,	273 }, /* 684    	*/
        {SubFunction , false, None  , 1, false, 0,	274 }, /* 685    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 686    	*/
        {SubFunction , false, None  , 1, false, 0,	275 }, /* 687    	*/
        {SubFunction , false, None  , 1, false, 0,	264 }, /* 688    	*/
        {SubFunction , false, None  , 1, false, 0,	266 }, /* 689    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 690    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 691    	*/
        {SubFunction , false, None  , 1, false, 0,	276 }, /* 692    	*/
        {SubFunction , false, None  , 1, false, 0,	265 }, /* 693    	*/
        {SubFunction , false, None  , 1, false, 0,	277 }, /* 694    	*/
        {SubFunction , false, None  , 1, false, 0,	276 }, /* 695    	*/
        {SubFunction , false, None  , 1, false, 0,	264 }, /* 696    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 697    	*/
        {SubFunction , false, None  , 1, false, 0,	16 }, /* 698    	*/
        {SubFunction , false, None  , 1, false, 0,	278 }, /* 699    	*/
        {SubFunction , false, None  , 1, false, 0,	276 }, /* 700    	*/
        {SubFunction , false, None  , 1, false, 0,	267 }, /* 701    	*/
        {SubFunction , false, None  , 1, false, 0,	279 }, /* 702    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 703    	*/
        {SubFunction , false, None  , 1, false, 0,	268 }, /* 704    	*/
        {SubFunction , false, None  , 1, false, 0,	280 }, /* 705    	*/
        {SubFunction , false, None  , 1, false, 0,	269 }, /* 706    	*/
        {SubFunction , false, None  , 1, false, 0,	274 }, /* 707    	*/
        {NumConstant , false, None  , 1, false, 0,	9 }, /* 708    	*/
        {SubFunction , false, None  , 1, false, 0,	281 }, /* 709    	*/
        {SubFunction , false, None  , 1, false, 0,	267 }, /* 710    	*/
        {SubFunction , false, None  , 1, false, 0,	265 }, /* 711    	*/
        {SubFunction , false, None  , 1, false, 0,	282 }, /* 712    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 713    	*/
        {SubFunction , false, None  , 1, false, 0,	268 }, /* 714    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 715    	*/
        {SubFunction , false, None  , 1, false, 0,	269 }, /* 716    	*/
        {SubFunction , false, None  , 1, false, 0,	283 }, /* 717    	*/
        {SubFunction , false, None  , 1, false, 0,	284 }, /* 718    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 719    	*/
        {ImmedHolder , false, None  , 1, false, 0,	0 }, /* 720    	*/
        {SubFunction , false, None  , 1, false, 0,	271 }, /* 721    	*/
        {cMul        , false, None  , 1, false, 2,	719 }, /* 722    	*/
        {SubFunction , false, None  , 1, false, 0,	285 }, /* 723    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 724    	*/
        {SubFunction , false, None  , 1, false, 0,	269 }, /* 725    	*/
        {SubFunction , false, None  , 1, false, 0,	283 }, /* 726    	*/
        {SubFunction , false, None  , 1, false, 0,	286 }, /* 727    	*/
        {SubFunction , false, None  , 1, false, 0,	274 }, /* 728    	*/
        {cMul        , false, None  , 1, false, 2,	719 }, /* 729    	*/
        {SubFunction , false, None  , 1, false, 0,	287 }, /* 730    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 731    	*/
        {SubFunction , false, None  , 1, false, 0,	288 }, /* 732    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 733    	*/
        {SubFunction , false, None  , 1, false, 0,	289 }, /* 734    	*/
        {SubFunction , false, None  , 1, false, 0,	290 }, /* 735    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 736    	*/
        {SubFunction , false, None  , 1, false, 0,	291 }, /* 737    	*/
        {SubFunction , false, None  , 1, false, 0,	292 }, /* 738    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 739    	*/
        {SubFunction , false, None  , 1, false, 0,	293 }, /* 740    	*/
        {SubFunction , false, None  , 1, false, 0,	294 }, /* 741    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 742    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 743    	*/
        {SubFunction , false, None  , 1, false, 0,	295 }, /* 744    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 745    	*/
        {SubFunction , false, None  , 1, false, 0,	296 }, /* 746    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 747    	*/
        {SubFunction , false, None  , 1, false, 0,	297 }, /* 748    	*/
        {SubFunction , false, None  , 1, false, 0,	298 }, /* 749    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 750    	*/
        {SubFunction , false, None  , 1, false, 0,	299 }, /* 751    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 752    	*/
        {SubFunction , false, None  , 1, false, 0,	300 }, /* 753    	*/
        {SubFunction , false, None  , 1, false, 0,	297 }, /* 754    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 755    	*/
        {SubFunction , false, None  , 1, false, 0,	300 }, /* 756    	*/
        {SubFunction , false, None  , 1, false, 0,	264 }, /* 757    	*/
        {SubFunction , false, None  , 1, false, 0,	301 }, /* 758    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 759    	*/
        {SubFunction , false, None  , 1, false, 0,	302 }, /* 760    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 761    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 762 "x"	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 763    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 764    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	1 }, /* 765    	*/
        {SubFunction , false, None  , 1, false, 0,	304 }, /* 766    	*/
        {cAdd        , false, None  , 1, false, 2,	764 }, /* 767    	*/
        {SubFunction , false, None  , 1, false, 0,	303 }, /* 768    	*/
        {SubFunction , false, None  , 1, false, 0,	305 }, /* 769    	*/
        {SubFunction , false, None  , 1, false, 0,	306 }, /* 770    	*/
        {ImmedHolder , false, None  , 1, false, 0,	1 }, /* 771    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 772    	*/
        {SubFunction , false, None  , 1, false, 0,	307 }, /* 773    	*/
        {NumConstant , false, None  , 1, false, 0,	0 }, /* 774    	*/
        {cLog        , false, Invert, 1, false, 1,	82 }, /* 775    	*/
        {SubFunction , false, None  , 1, false, 0,	308 }, /* 776    	*/
        {cMul        , false, None  , 1, false, 2,	774 }, /* 777    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 778    	*/
        {SubFunction , false, None  , 1, false, 0,	309 }, /* 779    	*/
        {SubFunction , false, None  , 1, false, 0,	310 }, /* 780    	*/
        {SubFunction , false, None  , 1, false, 0,	311 }, /* 781    	*/
        {SubFunction , false, None  , 1, false, 0,	312 }, /* 782    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 783    	*/
        {SubFunction , false, None  , 1, false, 0,	288 }, /* 784    	*/
        {SubFunction , false, None  , 1, false, 0,	313 }, /* 785    	*/
        {SubFunction , false, None  , 1, false, 0,	314 }, /* 786    	*/
        {SubFunction , false, None  , 1, false, 0,	315 }, /* 787    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 788    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 789    	*/
        {SubFunction , false, None  , 1, false, 0,	264 }, /* 790    	*/
        {SubFunction , false, None  , 1, false, 0,	316 }, /* 791    	*/
        {SubFunction , false, None  , 1, false, 0,	317 }, /* 792    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 793    	*/
        {SubFunction , false, None  , 1, false, 0,	289 }, /* 794    	*/
        {SubFunction , false, None  , 1, false, 0,	318 }, /* 795    	*/
        {SubFunction , false, None  , 1, false, 0,	319 }, /* 796    	*/
        {SubFunction , false, None  , 1, false, 0,	296 }, /* 797    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 798    	*/
        {SubFunction , false, None  , 1, false, 0,	296 }, /* 799    	*/
        {SubFunction , false, None  , 1, false, 0,	295 }, /* 800    	*/
        {NumConstant , false, None  , 1, false, 0,	6 }, /* 801    	*/
        {SubFunction , false, None  , 1, false, 0,	224 }, /* 802    	*/
        {SubFunction , false, None  , 1, false, 0,	320 }, /* 803    	*/
        {SubFunction , false, None  , 1, false, 0,	322 }, /* 804    	*/
        {SubFunction , false, None  , 1, false, 0,	321 }, /* 805    	*/
        {SubFunction , false, None  , 1, false, 0,	323 }, /* 806    	*/
        {SubFunction , false, None  , 1, false, 0,	299 }, /* 807    	*/
        {SubFunction , false, None  , 1, false, 0,	263 }, /* 808    	*/
        {SubFunction , false, None  , 1, false, 0,	299 }, /* 809    	*/
        {SubFunction , false, None  , 1, false, 0,	298 }, /* 810    	*/
        {NumConstant , false, None  , 1, false, 0,	10 }, /* 811    	*/
        {NamedHolder , false, None  , 1, false, 0,	0 }, /* 812 "x"	*/
        {SubFunction , false, None  , 1, false, 0,	324 }, /* 813    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 814    	*/
        {SubFunction , false, None  , 1, false, 0,	324 }, /* 815    	*/
        {NumConstant , false, None  , 1, false, 0,	3 }, /* 816    	*/
        {SubFunction , false, None  , 1, false, 0,	326 }, /* 817    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 818    	*/
        {SubFunction , false, None  , 1, false, 0,	325 }, /* 819    	*/
        {SubFunction , false, None  , 1, false, 0,	327 }, /* 820    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 821    	*/
        {SubFunction , false, None  , 1, false, 0,	299 }, /* 822    	*/
        {SubFunction , false, None  , 1, false, 0,	325 }, /* 823    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 824    	*/
        {SubFunction , false, None  , 1, false, 0,	326 }, /* 825    	*/
        {SubFunction , false, None  , 1, false, 0,	328 }, /* 826    	*/
        {SubFunction , false, None  , 1, false, 0,	329 }, /* 827    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 828    	*/
        {SubFunction , false, None  , 1, false, 0,	330 }, /* 829    	*/
        {SubFunction , false, None  , 1, false, 0,	332 }, /* 830    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 831    	*/
        {SubFunction , false, None  , 1, false, 0,	331 }, /* 832    	*/
        {SubFunction , false, None  , 1, false, 0,	333 }, /* 833    	*/
        {SubFunction , false, None  , 1, false, 0,	331 }, /* 834    	*/
        {NumConstant , false, None  , 1, false, 0,	1 }, /* 835    	*/
        {SubFunction , false, None  , 1, false, 0,	332 }, /* 836    	*/
        {SubFunction , false, None  , 1, false, 0,	334 }, /* 837    	*/
        {SubFunction , false, None  , 1, false, 0,	23 }, /* 838    	*/
        {ImmedHolder , false, None  , 1, false, Sign_Negative,	0 }, /* 839    	*/
        {SubFunction , false, None  , 1, false, 0,	335 }, /* 840    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 841    	*/
        {SubFunction , false, None  , 1, false, 0,	22 }, /* 842    	*/
        {ImmedHolder , false, None  , 1, false, Sign_Negative,	0 }, /* 843    	*/
        {SubFunction , false, None  , 1, false, 0,	336 }, /* 844    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 845    	*/
        {SubFunction , false, None  , 1, false, 0,	296 }, /* 846    	*/
        {ImmedHolder , false, None  , 1, false, Sign_Negative,	0 }, /* 847    	*/
        {SubFunction , false, None  , 1, false, 0,	337 }, /* 848    	*/
        {ImmedHolder , false, Negate, 1, false, 0,	0 }, /* 849    	*/
        {NumConstant , false, None  , 1, false, 0,	11 }, /* 850    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 851    	*/
        {SubFunction , false, None  , 1, false, 0,	338 }, /* 852    	*/
        {NumConstant , false, None  , 1, false, 0,	12 }, /* 853    	*/
        {RestHolder  , false, None  , 1, false, 0,	1 }, /* 854    	*/
        {SubFunction , false, None  , 1, false, 0,	339 }, /* 855    	*/
    };

    const MatchedParams mlist[] =
    {
        {PositionalParams, BalanceDontCare, 1, 0 }, /* 0 */
        {PositionalParams, BalanceDontCare, 2, 1 }, /* 1 */
        {PositionalParams, BalanceDontCare, 2, 3 }, /* 2 */
        {PositionalParams, BalanceDontCare, 2, 5 }, /* 3 */
        {PositionalParams, BalanceDontCare, 2, 7 }, /* 4 */
        {PositionalParams, BalanceDontCare, 1, 9 }, /* 5 */
        {PositionalParams, BalanceDontCare, 1, 10 }, /* 6 */
        {PositionalParams, BalanceDontCare, 2, 11 }, /* 7 */
        {PositionalParams, BalanceDontCare, 2, 13 }, /* 8 */
        {PositionalParams, BalanceDontCare, 2, 15 }, /* 9 */
        {PositionalParams, BalanceDontCare, 1, 17 }, /* 10 */
        {PositionalParams, BalanceDontCare, 1, 18 }, /* 11 */
        {PositionalParams, BalanceDontCare, 2, 19 }, /* 12 */
        {PositionalParams, BalanceDontCare, 2, 21 }, /* 13 */
        {PositionalParams, BalanceDontCare, 2, 23 }, /* 14 */
        {PositionalParams, BalanceDontCare, 2, 25 }, /* 15 */
        {PositionalParams, BalanceDontCare, 2, 27 }, /* 16 */
        {PositionalParams, BalanceDontCare, 1, 29 }, /* 17 */
        {PositionalParams, BalanceDontCare, 2, 30 }, /* 18 */
        {PositionalParams, BalanceDontCare, 1, 32 }, /* 19 */
        {PositionalParams, BalanceDontCare, 2, 33 }, /* 20 */
        {PositionalParams, BalanceDontCare, 2, 35 }, /* 21 */
        {PositionalParams, BalanceDontCare, 2, 37 }, /* 22 */
        {PositionalParams, BalanceDontCare, 3, 39 }, /* 23 */
        {PositionalParams, BalanceDontCare, 1, 42 }, /* 24 */
        {PositionalParams, BalanceDontCare, 2, 43 }, /* 25 */
        {PositionalParams, BalanceDontCare, 3, 45 }, /* 26 */
        {PositionalParams, BalanceDontCare, 1, 48 }, /* 27 */
        {PositionalParams, BalanceDontCare, 2, 49 }, /* 28 */
        {PositionalParams, BalanceDontCare, 2, 51 }, /* 29 */
        {PositionalParams, BalanceDontCare, 1, 53 }, /* 30 */
        {PositionalParams, BalanceDontCare, 2, 54 }, /* 31 */
        {PositionalParams, BalanceDontCare, 2, 56 }, /* 32 */
        {PositionalParams, BalanceDontCare, 2, 58 }, /* 33 */
        {PositionalParams, BalanceDontCare, 2, 60 }, /* 34 */
        {PositionalParams, BalanceDontCare, 1, 62 }, /* 35 */
        {PositionalParams, BalanceDontCare, 1, 63 }, /* 36 */
        {PositionalParams, BalanceDontCare, 1, 64 }, /* 37 */
        {AnyParams       , BalanceDontCare, 2, 65 }, /* 38 */
        {PositionalParams, BalanceDontCare, 1, 67 }, /* 39 */
        {PositionalParams, BalanceDontCare, 2, 68 }, /* 40 */
        {PositionalParams, BalanceDontCare, 1, 70 }, /* 41 */
        {PositionalParams, BalanceDontCare, 1, 71 }, /* 42 */
        {AnyParams       , BalanceDontCare, 2, 72 }, /* 43 */
        {PositionalParams, BalanceDontCare, 1, 74 }, /* 44 */
        {PositionalParams, BalanceDontCare, 1, 66 }, /* 45 */
        {PositionalParams, BalanceDontCare, 1, 75 }, /* 46 */
        {PositionalParams, BalanceDontCare, 1, 76 }, /* 47 */
        {AnyParams       , BalanceDontCare, 2, 77 }, /* 48 */
        {SelectedParams  , BalanceDontCare, 2, 79 }, /* 49 */
        {PositionalParams, BalanceDontCare, 1, 81 }, /* 50 */
        {PositionalParams, BalanceDontCare, 2, 82 }, /* 51 */
        {PositionalParams, BalanceDontCare, 1, 84 }, /* 52 */
        {PositionalParams, BalanceDontCare, 1, 85 }, /* 53 */
        {PositionalParams, BalanceDontCare, 1, 80 }, /* 54 */
        {PositionalParams, BalanceDontCare, 3, 86 }, /* 55 */
        {AnyParams       , BalanceDontCare, 2, 89 }, /* 56 */
        {PositionalParams, BalanceDontCare, 3, 91 }, /* 57 */
        {PositionalParams, BalanceDontCare, 3, 94 }, /* 58 */
        {PositionalParams, BalanceDontCare, 2, 97 }, /* 59 */
        {PositionalParams, BalanceDontCare, 1, 99 }, /* 60 */
        {PositionalParams, BalanceDontCare, 3, 100 }, /* 61 */
        {PositionalParams, BalanceDontCare, 3, 103 }, /* 62 */
        {PositionalParams, BalanceDontCare, 2, 106 }, /* 63 */
        {PositionalParams, BalanceDontCare, 1, 108 }, /* 64 */
        {PositionalParams, BalanceDontCare, 3, 109 }, /* 65 */
        {PositionalParams, BalanceDontCare, 3, 112 }, /* 66 */
        {PositionalParams, BalanceDontCare, 2, 115 }, /* 67 */
        {PositionalParams, BalanceDontCare, 1, 117 }, /* 68 */
        {PositionalParams, BalanceDontCare, 3, 118 }, /* 69 */
        {PositionalParams, BalanceDontCare, 3, 121 }, /* 70 */
        {PositionalParams, BalanceDontCare, 2, 124 }, /* 71 */
        {PositionalParams, BalanceDontCare, 1, 126 }, /* 72 */
        {AnyParams       , BalanceDontCare, 2, 127 }, /* 73 */
        {PositionalParams, BalanceDontCare, 3, 129 }, /* 74 */
        {PositionalParams, BalanceDontCare, 1, 128 }, /* 75 */
        {PositionalParams, BalanceDontCare, 3, 132 }, /* 76 */
        {PositionalParams, BalanceDontCare, 2, 135 }, /* 77 */
        {PositionalParams, BalanceDontCare, 1, 137 }, /* 78 */
        {PositionalParams, BalanceDontCare, 3, 138 }, /* 79 */
        {PositionalParams, BalanceDontCare, 3, 141 }, /* 80 */
        {PositionalParams, BalanceDontCare, 2, 144 }, /* 81 */
        {PositionalParams, BalanceDontCare, 1, 146 }, /* 82 */
        {PositionalParams, BalanceDontCare, 3, 147 }, /* 83 */
        {PositionalParams, BalanceDontCare, 3, 150 }, /* 84 */
        {PositionalParams, BalanceDontCare, 2, 153 }, /* 85 */
        {PositionalParams, BalanceDontCare, 1, 155 }, /* 86 */
        {PositionalParams, BalanceDontCare, 3, 156 }, /* 87 */
        {PositionalParams, BalanceDontCare, 3, 159 }, /* 88 */
        {PositionalParams, BalanceDontCare, 2, 162 }, /* 89 */
        {PositionalParams, BalanceDontCare, 1, 164 }, /* 90 */
        {AnyParams       , BalanceDontCare, 1, 0 }, /* 91 */
        {PositionalParams, BalanceDontCare, 3, 165 }, /* 92 */
        {PositionalParams, BalanceDontCare, 3, 168 }, /* 93 */
        {PositionalParams, BalanceDontCare, 1, 171 }, /* 94 */
        {PositionalParams, BalanceDontCare, 3, 172 }, /* 95 */
        {PositionalParams, BalanceDontCare, 3, 175 }, /* 96 */
        {PositionalParams, BalanceDontCare, 1, 178 }, /* 97 */
        {PositionalParams, BalanceDontCare, 2, 179 }, /* 98 */
        {PositionalParams, BalanceDontCare, 1, 181 }, /* 99 */
        {PositionalParams, BalanceDontCare, 1, 68 }, /* 100 */
        {PositionalParams, BalanceDontCare, 2, 182 }, /* 101 */
        {PositionalParams, BalanceDontCare, 1, 184 }, /* 102 */
        {PositionalParams, BalanceDontCare, 2, 185 }, /* 103 */
        {PositionalParams, BalanceDontCare, 1, 187 }, /* 104 */
        {PositionalParams, BalanceDontCare, 2, 188 }, /* 105 */
        {PositionalParams, BalanceDontCare, 1, 190 }, /* 106 */
        {PositionalParams, BalanceDontCare, 2, 191 }, /* 107 */
        {AnyParams       , BalanceDontCare, 2, 193 }, /* 108 */
        {PositionalParams, BalanceDontCare, 1, 195 }, /* 109 */
        {PositionalParams, BalanceDontCare, 1, 105 }, /* 110 */
        {PositionalParams, BalanceDontCare, 2, 196 }, /* 111 */
        {PositionalParams, BalanceDontCare, 2, 198 }, /* 112 */
        {PositionalParams, BalanceDontCare, 2, 200 }, /* 113 */
        {PositionalParams, BalanceDontCare, 1, 202 }, /* 114 */
        {AnyParams       , BalanceDontCare, 2, 203 }, /* 115 */
        {PositionalParams, BalanceDontCare, 1, 205 }, /* 116 */
        {PositionalParams, BalanceDontCare, 2, 206 }, /* 117 */
        {PositionalParams, BalanceDontCare, 1, 208 }, /* 118 */
        {AnyParams       , BalanceDontCare, 2, 0 }, /* 119 */
        {PositionalParams, BalanceDontCare, 2, 209 }, /* 120 */
        {AnyParams       , BalanceDontCare, 2, 211 }, /* 121 */
        {PositionalParams, BalanceDontCare, 2, 213 }, /* 122 */
        {PositionalParams, BalanceDontCare, 2, 215 }, /* 123 */
        {PositionalParams, BalanceDontCare, 2, 217 }, /* 124 */
        {PositionalParams, BalanceDontCare, 2, 219 }, /* 125 */
        {AnyParams       , BalanceDontCare, 2, 221 }, /* 126 */
        {PositionalParams, BalanceDontCare, 2, 223 }, /* 127 */
        {PositionalParams, BalanceDontCare, 2, 227 }, /* 128 */
        {AnyParams       , BalanceDontCare, 3, 229 }, /* 129 */
        {PositionalParams, BalanceDontCare, 2, 232 }, /* 130 */
        {SelectedParams  , BalanceDontCare, 2, 234 }, /* 131 */
        {PositionalParams, BalanceDontCare, 2, 236 }, /* 132 */
        {PositionalParams, BalanceDontCare, 2, 238 }, /* 133 */
        {PositionalParams, BalanceDontCare, 2, 240 }, /* 134 */
        {PositionalParams, BalanceDontCare, 2, 22 }, /* 135 */
        {PositionalParams, BalanceDontCare, 1, 167 }, /* 136 */
        {PositionalParams, BalanceDontCare, 2, 242 }, /* 137 */
        {AnyParams       , BalanceDontCare, 3, 244 }, /* 138 */
        {PositionalParams, BalanceDontCare, 2, 247 }, /* 139 */
        {PositionalParams, BalanceDontCare, 2, 249 }, /* 140 */
        {PositionalParams, BalanceDontCare, 2, 175 }, /* 141 */
        {PositionalParams, BalanceDontCare, 2, 251 }, /* 142 */
        {PositionalParams, BalanceDontCare, 2, 253 }, /* 143 */
        {PositionalParams, BalanceDontCare, 2, 166 }, /* 144 */
        {PositionalParams, BalanceDontCare, 2, 255 }, /* 145 */
        {PositionalParams, BalanceDontCare, 2, 257 }, /* 146 */
        {PositionalParams, BalanceDontCare, 2, 259 }, /* 147 */
        {PositionalParams, BalanceDontCare, 2, 261 }, /* 148 */
        {AnyParams       , BalanceDontCare, 2, 263 }, /* 149 */
        {PositionalParams, BalanceDontCare, 2, 265 }, /* 150 */
        {PositionalParams, BalanceDontCare, 2, 267 }, /* 151 */
        {PositionalParams, BalanceDontCare, 2, 269 }, /* 152 */
        {PositionalParams, BalanceDontCare, 1, 271 }, /* 153 */
        {PositionalParams, BalanceDontCare, 2, 272 }, /* 154 */
        {PositionalParams, BalanceDontCare, 2, 274 }, /* 155 */
        {PositionalParams, BalanceDontCare, 2, 276 }, /* 156 */
        {PositionalParams, BalanceDontCare, 1, 278 }, /* 157 */
        {PositionalParams, BalanceDontCare, 2, 279 }, /* 158 */
        {PositionalParams, BalanceDontCare, 1, 281 }, /* 159 */
        {PositionalParams, BalanceDontCare, 1, 282 }, /* 160 */
        {PositionalParams, BalanceDontCare, 2, 283 }, /* 161 */
        {PositionalParams, BalanceDontCare, 1, 285 }, /* 162 */
        {PositionalParams, BalanceDontCare, 2, 286 }, /* 163 */
        {SelectedParams  , BalanceDontCare, 2, 288 }, /* 164 */
        {AnyParams       , BalanceDontCare, 2, 290 }, /* 165 */
        {PositionalParams, BalanceDontCare, 2, 292 }, /* 166 */
        {PositionalParams, BalanceDontCare, 1, 294 }, /* 167 */
        {SelectedParams  , BalanceDontCare, 2, 295 }, /* 168 */
        {AnyParams       , BalanceDontCare, 2, 297 }, /* 169 */
        {PositionalParams, BalanceDontCare, 1, 289 }, /* 170 */
        {PositionalParams, BalanceDontCare, 1, 86 }, /* 171 */
        {AnyParams       , BalanceDontCare, 2, 299 }, /* 172 */
        {PositionalParams, BalanceDontCare, 1, 301 }, /* 173 */
        {PositionalParams, BalanceDontCare, 1, 302 }, /* 174 */
        {AnyParams       , BalanceDontCare, 2, 303 }, /* 175 */
        {PositionalParams, BalanceDontCare, 1, 12 }, /* 176 */
        {AnyParams       , BalanceDontCare, 2, 305 }, /* 177 */
        {AnyParams       , BalanceDontCare, 2, 307 }, /* 178 */
        {PositionalParams, BalanceDontCare, 1, 306 }, /* 179 */
        {PositionalParams, BalanceDontCare, 2, 309 }, /* 180 */
        {PositionalParams, BalanceDontCare, 2, 311 }, /* 181 */
        {PositionalParams, BalanceDontCare, 1, 313 }, /* 182 */
        {PositionalParams, BalanceDontCare, 2, 314 }, /* 183 */
        {AnyParams       , BalanceDontCare, 2, 316 }, /* 184 */
        {AnyParams       , BalanceDontCare, 2, 318 }, /* 185 */
        {PositionalParams, BalanceDontCare, 2, 322 }, /* 186 */
        {PositionalParams, BalanceDontCare, 2, 324 }, /* 187 */
        {PositionalParams, BalanceDontCare, 2, 326 }, /* 188 */
        {PositionalParams, BalanceDontCare, 2, 328 }, /* 189 */
        {PositionalParams, BalanceDontCare, 2, 330 }, /* 190 */
        {PositionalParams, BalanceDontCare, 2, 332 }, /* 191 */
        {PositionalParams, BalanceDontCare, 2, 334 }, /* 192 */
        {PositionalParams, BalanceDontCare, 2, 336 }, /* 193 */
        {PositionalParams, BalanceDontCare, 2, 338 }, /* 194 */
        {PositionalParams, BalanceDontCare, 1, 340 }, /* 195 */
        {PositionalParams, BalanceDontCare, 2, 341 }, /* 196 */
        {AnyParams       , BalanceDontCare, 2, 343 }, /* 197 */
        {AnyParams       , BalanceDontCare, 2, 345 }, /* 198 */
        {PositionalParams, BalanceDontCare, 2, 347 }, /* 199 */
        {PositionalParams, BalanceDontCare, 2, 349 }, /* 200 */
        {PositionalParams, BalanceDontCare, 2, 351 }, /* 201 */
        {PositionalParams, BalanceDontCare, 2, 353 }, /* 202 */
        {PositionalParams, BalanceDontCare, 2, 355 }, /* 203 */
        {PositionalParams, BalanceDontCare, 2, 357 }, /* 204 */
        {PositionalParams, BalanceDontCare, 2, 359 }, /* 205 */
        {PositionalParams, BalanceDontCare, 2, 361 }, /* 206 */
        {PositionalParams, BalanceDontCare, 2, 363 }, /* 207 */
        {PositionalParams, BalanceDontCare, 1, 365 }, /* 208 */
        {SelectedParams  , BalanceDontCare, 2, 366 }, /* 209 */
        {SelectedParams  , BalanceDontCare, 2, 368 }, /* 210 */
        {AnyParams       , BalanceDontCare, 2, 370 }, /* 211 */
        {PositionalParams, BalanceDontCare, 1, 372 }, /* 212 */
        {PositionalParams, BalanceDontCare, 1, 373 }, /* 213 */
        {SelectedParams  , BalanceDontCare, 3, 374 }, /* 214 */
        {AnyParams       , BalanceDontCare, 2, 377 }, /* 215 */
        {PositionalParams, BalanceDontCare, 2, 379 }, /* 216 */
        {PositionalParams, BalanceDontCare, 2, 381 }, /* 217 */
        {PositionalParams, BalanceDontCare, 1, 383 }, /* 218 */
        {PositionalParams, BalanceDontCare, 1, 384 }, /* 219 */
        {SelectedParams  , BalanceDontCare, 2, 385 }, /* 220 */
        {SelectedParams  , BalanceDontCare, 2, 387 }, /* 221 */
        {AnyParams       , BalanceDontCare, 2, 389 }, /* 222 */
        {PositionalParams, BalanceDontCare, 1, 391 }, /* 223 */
        {SelectedParams  , BalanceDontCare, 3, 392 }, /* 224 */
        {AnyParams       , BalanceDontCare, 2, 395 }, /* 225 */
        {PositionalParams, BalanceDontCare, 1, 397 }, /* 226 */
        {SelectedParams  , BalanceDontCare, 2, 398 }, /* 227 */
        {AnyParams       , BalanceDontCare, 2, 400 }, /* 228 */
        {PositionalParams, BalanceDontCare, 2, 404 }, /* 229 */
        {PositionalParams, BalanceDontCare, 2, 406 }, /* 230 */
        {PositionalParams, BalanceDontCare, 2, 408 }, /* 231 */
        {PositionalParams, BalanceDontCare, 1, 410 }, /* 232 */
        {SelectedParams  , BalanceDontCare, 3, 411 }, /* 233 */
        {AnyParams       , BalanceDontCare, 2, 414 }, /* 234 */
        {PositionalParams, BalanceDontCare, 2, 416 }, /* 235 */
        {PositionalParams, BalanceDontCare, 1, 418 }, /* 236 */
        {AnyParams       , BalanceDontCare, 2, 419 }, /* 237 */
        {PositionalParams, BalanceDontCare, 2, 421 }, /* 238 */
        {PositionalParams, BalanceDontCare, 1, 423 }, /* 239 */
        {SelectedParams  , BalanceDontCare, 3, 424 }, /* 240 */
        {AnyParams       , BalanceDontCare, 2, 427 }, /* 241 */
        {PositionalParams, BalanceDontCare, 2, 429 }, /* 242 */
        {PositionalParams, BalanceDontCare, 1, 431 }, /* 243 */
        {SelectedParams  , BalanceDontCare, 2, 432 }, /* 244 */
        {AnyParams       , BalanceDontCare, 2, 434 }, /* 245 */
        {AnyParams       , BalanceDontCare, 2, 436 }, /* 246 */
        {PositionalParams, BalanceDontCare, 2, 438 }, /* 247 */
        {PositionalParams, BalanceDontCare, 2, 440 }, /* 248 */
        {PositionalParams, BalanceDontCare, 1, 442 }, /* 249 */
        {AnyParams       , BalanceDontCare, 1, 443 }, /* 250 */
        {PositionalParams, BalanceDontCare, 2, 444 }, /* 251 */
        {PositionalParams, BalanceDontCare, 1, 446 }, /* 252 */
        {PositionalParams, BalanceDontCare, 2, 447 }, /* 253 */
        {PositionalParams, BalanceDontCare, 2, 449 }, /* 254 */
        {AnyParams       , BalanceDontCare, 3, 451 }, /* 255 */
        {AnyParams       , BalanceDontCare, 3, 454 }, /* 256 */
        {PositionalParams, BalanceDontCare, 2, 451 }, /* 257 */
        {PositionalParams, BalanceDontCare, 2, 457 }, /* 258 */
        {PositionalParams, BalanceDontCare, 2, 459 }, /* 259 */
        {PositionalParams, BalanceDontCare, 3, 461 }, /* 260 */
        {PositionalParams, BalanceDontCare, 2, 464 }, /* 261 */
        {SelectedParams  , BalanceDontCare, 2, 466 }, /* 262 */
        {SelectedParams  , BalanceDontCare, 2, 468 }, /* 263 */
        {SelectedParams  , BalanceDontCare, 2, 470 }, /* 264 */
        {PositionalParams, BalanceDontCare, 2, 234 }, /* 265 */
        {PositionalParams, BalanceDontCare, 1, 472 }, /* 266 */
        {AnyParams       , BalanceDontCare, 2, 473 }, /* 267 */
        {PositionalParams, BalanceDontCare, 2, 475 }, /* 268 */
        {AnyParams       , BalanceDontCare, 3, 477 }, /* 269 */
        {PositionalParams, BalanceDontCare, 1, 474 }, /* 270 */
        {PositionalParams, BalanceDontCare, 2, 480 }, /* 271 */
        {PositionalParams, BalanceDontCare, 2, 482 }, /* 272 */
        {PositionalParams, BalanceDontCare, 2, 168 }, /* 273 */
        {AnyParams       , BalanceDontCare, 2, 484 }, /* 274 */
        {PositionalParams, BalanceDontCare, 2, 486 }, /* 275 */
        {PositionalParams, BalanceDontCare, 1, 488 }, /* 276 */
        {AnyParams       , BalanceDontCare, 2, 489 }, /* 277 */
        {PositionalParams, BalanceDontCare, 2, 491 }, /* 278 */
        {PositionalParams, BalanceDontCare, 1, 493 }, /* 279 */
        {AnyParams       , BalanceDontCare, 2, 494 }, /* 280 */
        {AnyParams       , BalanceDontCare, 2, 496 }, /* 281 */
        {AnyParams       , BalanceDontCare, 2, 498 }, /* 282 */
        {PositionalParams, BalanceDontCare, 2, 500 }, /* 283 */
        {PositionalParams, BalanceDontCare, 2, 502 }, /* 284 */
        {PositionalParams, BalanceDontCare, 2, 504 }, /* 285 */
        {PositionalParams, BalanceDontCare, 1, 506 }, /* 286 */
        {AnyParams       , BalanceDontCare, 2, 507 }, /* 287 */
        {PositionalParams, BalanceDontCare, 2, 509 }, /* 288 */
        {PositionalParams, BalanceDontCare, 2, 511 }, /* 289 */
        {PositionalParams, BalanceDontCare, 1, 513 }, /* 290 */
        {AnyParams       , BalanceDontCare, 2, 514 }, /* 291 */
        {PositionalParams, BalanceDontCare, 2, 516 }, /* 292 */
        {PositionalParams, BalanceDontCare, 2, 518 }, /* 293 */
        {PositionalParams, BalanceDontCare, 1, 520 }, /* 294 */
        {PositionalParams, BalanceDontCare, 1, 521 }, /* 295 */
        {PositionalParams, BalanceDontCare, 2, 0 }, /* 296 */
        {PositionalParams, BalanceDontCare, 1, 95 }, /* 297 */
        {PositionalParams, BalanceDontCare, 1, 522 }, /* 298 */
        {PositionalParams, BalanceDontCare, 1, 523 }, /* 299 */
        {PositionalParams, BalanceDontCare, 1, 524 }, /* 300 */
        {PositionalParams, BalanceDontCare, 1, 525 }, /* 301 */
        {PositionalParams, BalanceDontCare, 1, 526 }, /* 302 */
        {PositionalParams, BalanceDontCare, 1, 527 }, /* 303 */
        {AnyParams       , BalanceMoreNeg , 2, 528 }, /* 304 */
        {PositionalParams, BalanceDontCare, 1, 530 }, /* 305 */
        {PositionalParams, BalanceDontCare, 2, 531 }, /* 306 */
        {PositionalParams, BalanceDontCare, 1, 533 }, /* 307 */
        {PositionalParams, BalanceDontCare, 1, 534 }, /* 308 */
        {PositionalParams, BalanceDontCare, 1, 535 }, /* 309 */
        {PositionalParams, BalanceDontCare, 2, 536 }, /* 310 */
        {PositionalParams, BalanceDontCare, 1, 538 }, /* 311 */
        {AnyParams       , BalanceDontCare, 1, 63 }, /* 312 */
        {PositionalParams, BalanceDontCare, 1, 539 }, /* 313 */
        {AnyParams       , BalanceDontCare, 1, 64 }, /* 314 */
        {AnyParams       , BalanceDontCare, 1, 540 }, /* 315 */
        {AnyParams       , BalanceDontCare, 1, 541 }, /* 316 */
        {AnyParams       , BalanceDontCare, 1, 542 }, /* 317 */
        {AnyParams       , BalanceDontCare, 1, 543 }, /* 318 */
        {AnyParams       , BalanceDontCare, 1, 544 }, /* 319 */
        {AnyParams       , BalanceDontCare, 1, 545 }, /* 320 */
        {AnyParams       , BalanceDontCare, 1, 546 }, /* 321 */
        {AnyParams       , BalanceDontCare, 1, 547 }, /* 322 */
        {AnyParams       , BalanceDontCare, 2, 548 }, /* 323 */
        {AnyParams       , BalanceDontCare, 2, 522 }, /* 324 */
        {AnyParams       , BalanceDontCare, 2, 550 }, /* 325 */
        {AnyParams       , BalanceDontCare, 2, 524 }, /* 326 */
        {AnyParams       , BalanceDontCare, 2, 552 }, /* 327 */
        {AnyParams       , BalanceDontCare, 2, 526 }, /* 328 */
        {AnyParams       , BalanceDontCare, 2, 554 }, /* 329 */
        {AnyParams       , BalanceDontCare, 2, 556 }, /* 330 */
        {AnyParams       , BalanceDontCare, 2, 558 }, /* 331 */
        {PositionalParams, BalanceDontCare, 2, 560 }, /* 332 */
        {PositionalParams, BalanceDontCare, 2, 562 }, /* 333 */
        {AnyParams       , BalanceDontCare, 3, 564 }, /* 334 */
        {PositionalParams, BalanceDontCare, 2, 564 }, /* 335 */
        {PositionalParams, BalanceDontCare, 1, 567 }, /* 336 */
        {AnyParams       , BalanceDontCare, 2, 568 }, /* 337 */
        {AnyParams       , BalanceDontCare, 2, 570 }, /* 338 */
        {AnyParams       , BalanceDontCare, 2, 572 }, /* 339 */
        {AnyParams       , BalanceDontCare, 2, 574 }, /* 340 */
        {PositionalParams, BalanceDontCare, 1, 576 }, /* 341 */
        {AnyParams       , BalanceDontCare, 2, 528 }, /* 342 */
        {PositionalParams, BalanceDontCare, 1, 577 }, /* 343 */
        {PositionalParams, BalanceDontCare, 2, 528 }, /* 344 */
        {PositionalParams, BalanceDontCare, 1, 578 }, /* 345 */
        {PositionalParams, BalanceDontCare, 1, 579 }, /* 346 */
        {PositionalParams, BalanceDontCare, 1, 580 }, /* 347 */
        {PositionalParams, BalanceDontCare, 2, 139 }, /* 348 */
        {PositionalParams, BalanceDontCare, 2, 142 }, /* 349 */
        {PositionalParams, BalanceDontCare, 1, 581 }, /* 350 */
        {PositionalParams, BalanceDontCare, 1, 16 }, /* 351 */
        {AnyParams       , BalanceDontCare, 2, 481 }, /* 352 */
        {PositionalParams, BalanceDontCare, 1, 582 }, /* 353 */
        {SelectedParams  , BalanceDontCare, 2, 3 }, /* 354 */
        {PositionalParams, BalanceDontCare, 2, 583 }, /* 355 */
        {SelectedParams  , BalanceDontCare, 2, 585 }, /* 356 */
        {PositionalParams, BalanceDontCare, 1, 587 }, /* 357 */
        {PositionalParams, BalanceDontCare, 1, 588 }, /* 358 */
        {SelectedParams  , BalanceDontCare, 2, 11 }, /* 359 */
        {PositionalParams, BalanceDontCare, 2, 589 }, /* 360 */
        {SelectedParams  , BalanceDontCare, 2, 591 }, /* 361 */
        {PositionalParams, BalanceDontCare, 1, 593 }, /* 362 */
        {SelectedParams  , BalanceDontCare, 2, 19 }, /* 363 */
        {SelectedParams  , BalanceDontCare, 2, 21 }, /* 364 */
        {SelectedParams  , BalanceDontCare, 2, 594 }, /* 365 */
        {PositionalParams, BalanceDontCare, 2, 596 }, /* 366 */
        {SelectedParams  , BalanceDontCare, 2, 598 }, /* 367 */
        {PositionalParams, BalanceDontCare, 1, 600 }, /* 368 */
        {PositionalParams, BalanceDontCare, 2, 601 }, /* 369 */
        {PositionalParams, BalanceDontCare, 1, 603 }, /* 370 */
        {SelectedParams  , BalanceDontCare, 2, 604 }, /* 371 */
        {SelectedParams  , BalanceDontCare, 2, 606 }, /* 372 */
        {SelectedParams  , BalanceDontCare, 2, 608 }, /* 373 */
        {PositionalParams, BalanceDontCare, 2, 610 }, /* 374 */
        {SelectedParams  , BalanceDontCare, 2, 612 }, /* 375 */
        {PositionalParams, BalanceDontCare, 1, 614 }, /* 376 */
        {PositionalParams, BalanceDontCare, 2, 604 }, /* 377 */
        {PositionalParams, BalanceDontCare, 1, 615 }, /* 378 */
        {PositionalParams, BalanceDontCare, 2, 616 }, /* 379 */
        {PositionalParams, BalanceDontCare, 1, 618 }, /* 380 */
        {PositionalParams, BalanceDontCare, 1, 8 }, /* 381 */
        {PositionalParams, BalanceDontCare, 2, 619 }, /* 382 */
        {PositionalParams, BalanceDontCare, 1, 621 }, /* 383 */
        {PositionalParams, BalanceDontCare, 1, 622 }, /* 384 */
        {PositionalParams, BalanceDontCare, 2, 623 }, /* 385 */
        {PositionalParams, BalanceDontCare, 1, 625 }, /* 386 */
        {PositionalParams, BalanceDontCare, 2, 626 }, /* 387 */
        {PositionalParams, BalanceDontCare, 1, 628 }, /* 388 */
        {PositionalParams, BalanceDontCare, 2, 629 }, /* 389 */
        {PositionalParams, BalanceDontCare, 1, 631 }, /* 390 */
        {PositionalParams, BalanceDontCare, 2, 632 }, /* 391 */
        {PositionalParams, BalanceDontCare, 1, 634 }, /* 392 */
        {SelectedParams  , BalanceDontCare, 2, 635 }, /* 393 */
        {PositionalParams, BalanceDontCare, 2, 639 }, /* 394 */
        {PositionalParams, BalanceDontCare, 1, 641 }, /* 395 */
        {PositionalParams, BalanceDontCare, 2, 644 }, /* 396 */
        {PositionalParams, BalanceDontCare, 3, 646 }, /* 397 */
        {PositionalParams, BalanceDontCare, 1, 649 }, /* 398 */
        {SelectedParams  , BalanceDontCare, 2, 650 }, /* 399 */
        {PositionalParams, BalanceDontCare, 3, 652 }, /* 400 */
        {PositionalParams, BalanceDontCare, 1, 655 }, /* 401 */
        {AnyParams       , BalanceDontCare, 2, 656 }, /* 402 */
        {PositionalParams, BalanceDontCare, 2, 658 }, /* 403 */
        {PositionalParams, BalanceDontCare, 1, 660 }, /* 404 */
        {SelectedParams  , BalanceDontCare, 2, 661 }, /* 405 */
        {AnyParams       , BalanceDontCare, 2, 663 }, /* 406 */
        {PositionalParams, BalanceDontCare, 1, 657 }, /* 407 */
        {SelectedParams  , BalanceDontCare, 2, 665 }, /* 408 */
        {AnyParams       , BalanceDontCare, 2, 667 }, /* 409 */
        {PositionalParams, BalanceDontCare, 1, 40 }, /* 410 */
        {PositionalParams, BalanceDontCare, 2, 669 }, /* 411 */
        {PositionalParams, BalanceDontCare, 2, 671 }, /* 412 */
        {AnyParams       , BalanceDontCare, 2, 673 }, /* 413 */
        {PositionalParams, BalanceDontCare, 2, 675 }, /* 414 */
        {PositionalParams, BalanceDontCare, 1, 677 }, /* 415 */
        {PositionalParams, BalanceDontCare, 2, 678 }, /* 416 */
        {PositionalParams, BalanceDontCare, 1, 680 }, /* 417 */
        {SelectedParams  , BalanceDontCare, 2, 681 }, /* 418 */
        {AnyParams       , BalanceDontCare, 2, 683 }, /* 419 */
        {PositionalParams, BalanceDontCare, 2, 685 }, /* 420 */
        {PositionalParams, BalanceDontCare, 1, 687 }, /* 421 */
        {AnyParams       , BalanceDontCare, 2, 688 }, /* 422 */
        {PositionalParams, BalanceDontCare, 1, 656 }, /* 423 */
        {SelectedParams  , BalanceDontCare, 2, 690 }, /* 424 */
        {AnyParams       , BalanceDontCare, 2, 692 }, /* 425 */
        {PositionalParams, BalanceDontCare, 2, 665 }, /* 426 */
        {PositionalParams, BalanceDontCare, 1, 694 }, /* 427 */
        {AnyParams       , BalanceDontCare, 2, 695 }, /* 428 */
        {PositionalParams, BalanceDontCare, 2, 697 }, /* 429 */
        {PositionalParams, BalanceDontCare, 1, 699 }, /* 430 */
        {AnyParams       , BalanceDontCare, 2, 700 }, /* 431 */
        {PositionalParams, BalanceDontCare, 2, 661 }, /* 432 */
        {PositionalParams, BalanceDontCare, 1, 702 }, /* 433 */
        {SelectedParams  , BalanceDontCare, 2, 703 }, /* 434 */
        {AnyParams       , BalanceDontCare, 2, 705 }, /* 435 */
        {PositionalParams, BalanceDontCare, 2, 707 }, /* 436 */
        {PositionalParams, BalanceDontCare, 1, 709 }, /* 437 */
        {AnyParams       , BalanceDontCare, 2, 710 }, /* 438 */
        {PositionalParams, BalanceDontCare, 2, 690 }, /* 439 */
        {PositionalParams, BalanceDontCare, 1, 712 }, /* 440 */
        {SelectedParams  , BalanceDontCare, 2, 713 }, /* 441 */
        {SelectedParams  , BalanceDontCare, 2, 715 }, /* 442 */
        {AnyParams       , BalanceDontCare, 2, 717 }, /* 443 */
        {PositionalParams, BalanceDontCare, 2, 721 }, /* 444 */
        {PositionalParams, BalanceDontCare, 1, 723 }, /* 445 */
        {SelectedParams  , BalanceDontCare, 2, 724 }, /* 446 */
        {AnyParams       , BalanceDontCare, 2, 726 }, /* 447 */
        {PositionalParams, BalanceDontCare, 2, 728 }, /* 448 */
        {PositionalParams, BalanceDontCare, 1, 730 }, /* 449 */
        {PositionalParams, BalanceDontCare, 1, 265 }, /* 450 */
        {AnyParams       , BalanceDontCare, 2, 731 }, /* 451 */
        {AnyParams       , BalanceDontCare, 2, 733 }, /* 452 */
        {PositionalParams, BalanceDontCare, 1, 735 }, /* 453 */
        {AnyParams       , BalanceDontCare, 2, 736 }, /* 454 */
        {PositionalParams, BalanceDontCare, 1, 738 }, /* 455 */
        {AnyParams       , BalanceDontCare, 2, 739 }, /* 456 */
        {PositionalParams, BalanceDontCare, 1, 741 }, /* 457 */
        {PositionalParams, BalanceDontCare, 2, 742 }, /* 458 */
        {AnyParams       , BalanceDontCare, 2, 744 }, /* 459 */
        {PositionalParams, BalanceDontCare, 2, 746 }, /* 460 */
        {PositionalParams, BalanceDontCare, 1, 748 }, /* 461 */
        {PositionalParams, BalanceDontCare, 2, 696 }, /* 462 */
        {AnyParams       , BalanceDontCare, 2, 749 }, /* 463 */
        {PositionalParams, BalanceDontCare, 2, 751 }, /* 464 */
        {PositionalParams, BalanceDontCare, 1, 753 }, /* 465 */
        {AnyParams       , BalanceDontCare, 2, 754 }, /* 466 */
        {PositionalParams, BalanceDontCare, 1, 49 }, /* 467 */
        {AnyParams       , BalanceDontCare, 2, 756 }, /* 468 */
        {SelectedParams  , BalanceDontCare, 2, 314 }, /* 469 */
        {SelectedParams  , BalanceDontCare, 2, 758 }, /* 470 */
        {PositionalParams, BalanceDontCare, 2, 760 }, /* 471 */
        {SelectedParams  , BalanceDontCare, 2, 762 }, /* 472 */
        {SelectedParams  , BalanceDontCare, 2, 766 }, /* 473 */
        {AnyParams       , BalanceDontCare, 2, 768 }, /* 474 */
        {PositionalParams, BalanceDontCare, 2, 432 }, /* 475 */
        {PositionalParams, BalanceDontCare, 3, 770 }, /* 476 */
        {PositionalParams, BalanceDontCare, 1, 773 }, /* 477 */
        {PositionalParams, BalanceDontCare, 2, 776 }, /* 478 */
        {PositionalParams, BalanceDontCare, 2, 778 }, /* 479 */
        {PositionalParams, BalanceDontCare, 1, 780 }, /* 480 */
        {AnyParams       , BalanceDontCare, 2, 51 }, /* 481 */
        {PositionalParams, BalanceDontCare, 1, 746 }, /* 482 */
        {AnyParams       , BalanceDontCare, 2, 82 }, /* 483 */
        {PositionalParams, BalanceDontCare, 1, 781 }, /* 484 */
        {PositionalParams, BalanceDontCare, 2, 782 }, /* 485 */
        {AnyParams       , BalanceDontCare, 2, 784 }, /* 486 */
        {PositionalParams, BalanceDontCare, 2, 263 }, /* 487 */
        {PositionalParams, BalanceDontCare, 1, 786 }, /* 488 */
        {PositionalParams, BalanceDontCare, 1, 787 }, /* 489 */
        {PositionalParams, BalanceDontCare, 2, 788 }, /* 490 */
        {AnyParams       , BalanceDontCare, 2, 790 }, /* 491 */
        {PositionalParams, BalanceDontCare, 1, 751 }, /* 492 */
        {PositionalParams, BalanceDontCare, 2, 792 }, /* 493 */
        {AnyParams       , BalanceDontCare, 2, 794 }, /* 494 */
        {PositionalParams, BalanceDontCare, 1, 796 }, /* 495 */
        {AnyParams       , BalanceDontCare, 2, 797 }, /* 496 */
        {PositionalParams, BalanceDontCare, 1, 51 }, /* 497 */
        {AnyParams       , BalanceDontCare, 2, 799 }, /* 498 */
        {PositionalParams, BalanceDontCare, 1, 52 }, /* 499 */
        {SelectedParams  , BalanceDontCare, 2, 801 }, /* 500 */
        {PositionalParams, BalanceDontCare, 1, 803 }, /* 501 */
        {AnyParams       , BalanceDontCare, 1, 66 }, /* 502 */
        {PositionalParams, BalanceDontCare, 1, 804 }, /* 503 */
        {AnyParams       , BalanceDontCare, 2, 805 }, /* 504 */
        {AnyParams       , BalanceDontCare, 2, 807 }, /* 505 */
        {AnyParams       , BalanceDontCare, 2, 809 }, /* 506 */
        {PositionalParams, BalanceDontCare, 1, 791 }, /* 507 */
        {PositionalParams, BalanceDontCare, 2, 811 }, /* 508 */
        {SelectedParams  , BalanceDontCare, 2, 813 }, /* 509 */
        {SelectedParams  , BalanceDontCare, 2, 815 }, /* 510 */
        {PositionalParams, BalanceDontCare, 2, 817 }, /* 511 */
        {AnyParams       , BalanceDontCare, 2, 819 }, /* 512 */
        {PositionalParams, BalanceDontCare, 2, 821 }, /* 513 */
        {PositionalParams, BalanceDontCare, 2, 823 }, /* 514 */
        {AnyParams       , BalanceDontCare, 2, 825 }, /* 515 */
        {PositionalParams, BalanceDontCare, 2, 827 }, /* 516 */
        {PositionalParams, BalanceDontCare, 1, 829 }, /* 517 */
        {SelectedParams  , BalanceDontCare, 2, 54 }, /* 518 */
        {SelectedParams  , BalanceDontCare, 2, 56 }, /* 519 */
        {PositionalParams, BalanceDontCare, 2, 830 }, /* 520 */
        {AnyParams       , BalanceDontCare, 2, 832 }, /* 521 */
        {PositionalParams, BalanceDontCare, 2, 834 }, /* 522 */
        {AnyParams       , BalanceDontCare, 2, 836 }, /* 523 */
        {PositionalParams, BalanceDontCare, 2, 838 }, /* 524 */
        {PositionalParams, BalanceDontCare, 2, 840 }, /* 525 */
        {PositionalParams, BalanceDontCare, 2, 842 }, /* 526 */
        {PositionalParams, BalanceDontCare, 2, 844 }, /* 527 */
        {PositionalParams, BalanceDontCare, 2, 846 }, /* 528 */
        {PositionalParams, BalanceDontCare, 2, 848 }, /* 529 */
        {AnyParams       , BalanceDontCare, 2, 850 }, /* 530 */
        {PositionalParams, BalanceDontCare, 1, 852 }, /* 531 */
        {AnyParams       , BalanceDontCare, 2, 853 }, /* 532 */
        {PositionalParams, BalanceDontCare, 1, 855 }, /* 533 */
    };

    const Function flist[] =
    {
        {cPow        , 1 }, /* 0 */
        {cAdd        , 2 }, /* 1 */
        {cPow        , 3 }, /* 2 */
        {cAdd        , 4 }, /* 3 */
        {cLog        , 5 }, /* 4 */
        {cAdd        , 7 }, /* 5 */
        {cPow        , 8 }, /* 6 */
        {cAdd        , 9 }, /* 7 */
        {cLog        , 10 }, /* 8 */
        {cAdd        , 12 }, /* 9 */
        {cMul        , 13 }, /* 10 */
        {cAdd        , 14 }, /* 11 */
        {cPow        , 15 }, /* 12 */
        {cMul        , 16 }, /* 13 */
        {cLog        , 17 }, /* 14 */
        {cMul        , 18 }, /* 15 */
        {cPow        , 20 }, /* 16 */
        {cPow        , 21 }, /* 17 */
        {cAdd        , 22 }, /* 18 */
        {cMul        , 23 }, /* 19 */
        {cAdd        , 25 }, /* 20 */
        {cMul        , 26 }, /* 21 */
        {cSin        , 0 }, /* 22 */
        {cCos        , 0 }, /* 23 */
        {cPow        , 28 }, /* 24 */
        {cMul        , 29 }, /* 25 */
        {cAdd        , 31 }, /* 26 */
        {cAdd        , 32 }, /* 27 */
        {cPow        , 33 }, /* 28 */
        {cMul        , 34 }, /* 29 */
        {cNot        , 0 }, /* 30 */
        {cNotNot     , 0 }, /* 31 */
        {cMul        , 38 }, /* 32 */
        {cAbs        , 0 }, /* 33 */
        {cMul        , 40 }, /* 34 */
        {cAcos       , 0 }, /* 35 */
        {cAdd        , 43 }, /* 36 */
        {cAdd        , 45 }, /* 37 */
        {cSin        , 46 }, /* 38 */
        {cMul        , 48 }, /* 39 */
        {cAdd        , 49 }, /* 40 */
        {cMul        , 51 }, /* 41 */
        {cSin        , 52 }, /* 42 */
        {cAdd        , 56 }, /* 43 */
        {cIf         , 58 }, /* 44 */
        {cAdd        , 59 }, /* 45 */
        {cMul        , 56 }, /* 46 */
        {cMul        , 45 }, /* 47 */
        {cIf         , 62 }, /* 48 */
        {cMul        , 63 }, /* 49 */
        {cAnd        , 56 }, /* 50 */
        {cAnd        , 45 }, /* 51 */
        {cIf         , 66 }, /* 52 */
        {cAnd        , 67 }, /* 53 */
        {cOr         , 56 }, /* 54 */
        {cOr         , 45 }, /* 55 */
        {cIf         , 70 }, /* 56 */
        {cOr         , 71 }, /* 57 */
        {cAdd        , 73 }, /* 58 */
        {cAdd        , 75 }, /* 59 */
        {cIf         , 76 }, /* 60 */
        {cAdd        , 77 }, /* 61 */
        {cMul        , 73 }, /* 62 */
        {cMul        , 75 }, /* 63 */
        {cIf         , 80 }, /* 64 */
        {cMul        , 81 }, /* 65 */
        {cAnd        , 73 }, /* 66 */
        {cAnd        , 75 }, /* 67 */
        {cIf         , 84 }, /* 68 */
        {cAnd        , 85 }, /* 69 */
        {cOr         , 73 }, /* 70 */
        {cOr         , 75 }, /* 71 */
        {cIf         , 88 }, /* 72 */
        {cOr         , 89 }, /* 73 */
        {cNot        , 91 }, /* 74 */
        {cIf         , 93 }, /* 75 */
        {cNotNot     , 91 }, /* 76 */
        {cIf         , 96 }, /* 77 */
        {cPow        , 98 }, /* 78 */
        {cLog        , 100 }, /* 79 */
        {cMul        , 101 }, /* 80 */
        {cPow        , 103 }, /* 81 */
        {cLog        , 0 }, /* 82 */
        {cMul        , 105 }, /* 83 */
        {cPow        , 107 }, /* 84 */
        {cMul        , 108 }, /* 85 */
        {cLog        , 110 }, /* 86 */
        {cMul        , 111 }, /* 87 */
        {cAdd        , 112 }, /* 88 */
        {cMul        , 113 }, /* 89 */
        {cMul        , 115 }, /* 90 */
        {cAdd        , 117 }, /* 91 */
        {cMul        , 121 }, /* 92 */
        {cMul        , 126 }, /* 93 */
        {cMul        , 129 }, /* 94 */
        {cAdd        , 131 }, /* 95 */
        {cMul        , 133 }, /* 96 */
        {cLog        , 136 }, /* 97 */
        {cPow        , 137 }, /* 98 */
        {cMul        , 138 }, /* 99 */
        {cPow        , 142 }, /* 100 */
        {cMul        , 144 }, /* 101 */
        {cPow        , 146 }, /* 102 */
        {cMul        , 149 }, /* 103 */
        {cPow        , 151 }, /* 104 */
        {cMul        , 152 }, /* 105 */
        {cPow        , 154 }, /* 106 */
        {cAsin       , 0 }, /* 107 */
        {cCos        , 46 }, /* 108 */
        {cMul        , 158 }, /* 109 */
        {cCos        , 52 }, /* 110 */
        {cMul        , 161 }, /* 111 */
        {cPow        , 163 }, /* 112 */
        {cMul        , 164 }, /* 113 */
        {cPow        , 166 }, /* 114 */
        {cMul        , 168 }, /* 115 */
        {cLog        , 171 }, /* 116 */
        {cMul        , 141 }, /* 117 */
        {cLog        , 173 }, /* 118 */
        {cMul        , 177 }, /* 119 */
        {cMul        , 179 }, /* 120 */
        {cAdd        , 180 }, /* 121 */
        {cMul        , 181 }, /* 122 */
        {cPow        , 183 }, /* 123 */
        {cMul        , 184 }, /* 124 */
        {cPow        , 186 }, /* 125 */
        {cAdd        , 187 }, /* 126 */
        {cPow        , 188 }, /* 127 */
        {cMul        , 189 }, /* 128 */
        {cAdd        , 190 }, /* 129 */
        {cPow        , 191 }, /* 130 */
        {cMul        , 192 }, /* 131 */
        {cAdd        , 193 }, /* 132 */
        {cMul        , 194 }, /* 133 */
        {cPow        , 196 }, /* 134 */
        {cMul        , 197 }, /* 135 */
        {cPow        , 199 }, /* 136 */
        {cAdd        , 200 }, /* 137 */
        {cPow        , 201 }, /* 138 */
        {cMul        , 202 }, /* 139 */
        {cAdd        , 203 }, /* 140 */
        {cPow        , 204 }, /* 141 */
        {cMul        , 205 }, /* 142 */
        {cAdd        , 206 }, /* 143 */
        {cMul        , 207 }, /* 144 */
        {cCos        , 171 }, /* 145 */
        {cMul        , 209 }, /* 146 */
        {cSin        , 171 }, /* 147 */
        {cMul        , 210 }, /* 148 */
        {cAdd        , 141 }, /* 149 */
        {cCos        , 212 }, /* 150 */
        {cMul        , 214 }, /* 151 */
        {cMul        , 216 }, /* 152 */
        {cAdd        , 217 }, /* 153 */
        {cCos        , 218 }, /* 154 */
        {cMul        , 220 }, /* 155 */
        {cMul        , 221 }, /* 156 */
        {cSin        , 212 }, /* 157 */
        {cMul        , 224 }, /* 158 */
        {cSin        , 218 }, /* 159 */
        {cMul        , 227 }, /* 160 */
        {cMul        , 229 }, /* 161 */
        {cAdd        , 230 }, /* 162 */
        {cMul        , 231 }, /* 163 */
        {cMul        , 233 }, /* 164 */
        {cMul        , 235 }, /* 165 */
        {cMul        , 238 }, /* 166 */
        {cMul        , 240 }, /* 167 */
        {cMul        , 242 }, /* 168 */
        {cMul        , 244 }, /* 169 */
        {cAdd        , 247 }, /* 170 */
        {cMul        , 248 }, /* 171 */
        {cMul        , 251 }, /* 172 */
        {cPow        , 253 }, /* 173 */
        {cPow        , 254 }, /* 174 */
        {cMul        , 255 }, /* 175 */
        {cAdd        , 257 }, /* 176 */
        {cPow        , 258 }, /* 177 */
        {cAdd        , 259 }, /* 178 */
        {cMul        , 260 }, /* 179 */
        {cMul        , 262 }, /* 180 */
        {cAdd        , 263 }, /* 181 */
        {cAdd        , 265 }, /* 182 */
        {cMul        , 267 }, /* 183 */
        {cPow        , 268 }, /* 184 */
        {cMul        , 270 }, /* 185 */
        {cPow        , 271 }, /* 186 */
        {cPow        , 141 }, /* 187 */
        {cPow        , 273 }, /* 188 */
        {cAdd        , 144 }, /* 189 */
        {cPow        , 275 }, /* 190 */
        {cPow        , 144 }, /* 191 */
        {cPow        , 278 }, /* 192 */
        {cMul        , 280 }, /* 193 */
        {cAdd        , 281 }, /* 194 */
        {cMul        , 283 }, /* 195 */
        {cMul        , 284 }, /* 196 */
        {cAdd        , 285 }, /* 197 */
        {cAdd        , 280 }, /* 198 */
        {cMul        , 288 }, /* 199 */
        {cAdd        , 289 }, /* 200 */
        {cAdd        , 292 }, /* 201 */
        {cPow        , 293 }, /* 202 */
        {cPow        , 251 }, /* 203 */
        {cEqual      , 257 }, /* 204 */
        {cNEqual     , 257 }, /* 205 */
        {cLess       , 257 }, /* 206 */
        {cGreaterOrEq, 257 }, /* 207 */
        {cLessOrEq   , 257 }, /* 208 */
        {cGreater    , 257 }, /* 209 */
        {cAnd        , 304 }, /* 210 */
        {cOr         , 306 }, /* 211 */
        {cOr         , 304 }, /* 212 */
        {cAnd        , 306 }, /* 213 */
        {cOr         , 310 }, /* 214 */
        {cEqual      , 332 }, /* 215 */
        {cEqual      , 333 }, /* 216 */
        {cAnd        , 310 }, /* 217 */
        {cNot        , 37 }, /* 218 */
        {cAnd        , 342 }, /* 219 */
        {cAnd        , 344 }, /* 220 */
        {cOr         , 342 }, /* 221 */
        {cOr         , 344 }, /* 222 */
        {cAsinh      , 0 }, /* 223 */
        {cMul        , 352 }, /* 224 */
        {cAdd        , 354 }, /* 225 */
        {cPow        , 355 }, /* 226 */
        {cAdd        , 356 }, /* 227 */
        {cAcosh      , 0 }, /* 228 */
        {cAdd        , 359 }, /* 229 */
        {cPow        , 360 }, /* 230 */
        {cAdd        , 361 }, /* 231 */
        {cAdd        , 363 }, /* 232 */
        {cMul        , 364 }, /* 233 */
        {cAdd        , 365 }, /* 234 */
        {cPow        , 366 }, /* 235 */
        {cMul        , 367 }, /* 236 */
        {cAtanh      , 0 }, /* 237 */
        {cMul        , 369 }, /* 238 */
        {cMul        , 371 }, /* 239 */
        {cAdd        , 372 }, /* 240 */
        {cAdd        , 373 }, /* 241 */
        {cPow        , 374 }, /* 242 */
        {cMul        , 375 }, /* 243 */
        {cMul        , 377 }, /* 244 */
        {cAtanh      , 378 }, /* 245 */
        {cMul        , 379 }, /* 246 */
        {cSinh       , 110 }, /* 247 */
        {cMul        , 382 }, /* 248 */
        {cAtan       , 0 }, /* 249 */
        {cTan        , 110 }, /* 250 */
        {cMul        , 385 }, /* 251 */
        {cTanh       , 110 }, /* 252 */
        {cMul        , 387 }, /* 253 */
        {cPow        , 389 }, /* 254 */
        {cPow        , 391 }, /* 255 */
        {cPow        , 377 }, /* 256 */
        {cMul        , 394 }, /* 257 */
        {cSinh       , 395 }, /* 258 */
        {cPow        , 396 }, /* 259 */
        {cMul        , 397 }, /* 260 */
        {cCosh       , 395 }, /* 261 */
        {cMul        , 400 }, /* 262 */
        {cCosh       , 0 }, /* 263 */
        {cSinh       , 0 }, /* 264 */
        {cPow        , 403 }, /* 265 */
        {cMul        , 405 }, /* 266 */
        {cMul        , 408 }, /* 267 */
        {cPow        , 411 }, /* 268 */
        {cPow        , 412 }, /* 269 */
        {cMul        , 414 }, /* 270 */
        {cCosh       , 415 }, /* 271 */
        {cMul        , 416 }, /* 272 */
        {cMul        , 418 }, /* 273 */
        {cSinh       , 415 }, /* 274 */
        {cMul        , 420 }, /* 275 */
        {cMul        , 424 }, /* 276 */
        {cMul        , 426 }, /* 277 */
        {cMul        , 429 }, /* 278 */
        {cMul        , 432 }, /* 279 */
        {cMul        , 434 }, /* 280 */
        {cMul        , 436 }, /* 281 */
        {cMul        , 439 }, /* 282 */
        {cMul        , 441 }, /* 283 */
        {cMul        , 442 }, /* 284 */
        {cMul        , 444 }, /* 285 */
        {cMul        , 446 }, /* 286 */
        {cMul        , 448 }, /* 287 */
        {cSin        , 450 }, /* 288 */
        {cSinh       , 450 }, /* 289 */
        {cSinh       , 52 }, /* 290 */
        {cTan        , 450 }, /* 291 */
        {cTan        , 52 }, /* 292 */
        {cTanh       , 450 }, /* 293 */
        {cTanh       , 52 }, /* 294 */
        {cPow        , 458 }, /* 295 */
        {cTan        , 0 }, /* 296 */
        {cPow        , 460 }, /* 297 */
        {cPow        , 462 }, /* 298 */
        {cTanh       , 0 }, /* 299 */
        {cPow        , 464 }, /* 300 */
        {cMul        , 469 }, /* 301 */
        {cAdd        , 470 }, /* 302 */
        {cPow        , 471 }, /* 303 */
        {cMul        , 472 }, /* 304 */
        {cAdd        , 473 }, /* 305 */
        {cMul        , 475 }, /* 306 */
        {cAdd        , 476 }, /* 307 */
        {cAtanh      , 477 }, /* 308 */
        {cMul        , 478 }, /* 309 */
        {cPow        , 479 }, /* 310 */
        {cMul        , 483 }, /* 311 */
        {cCos        , 484 }, /* 312 */
        {cPow        , 485 }, /* 313 */
        {cMul        , 487 }, /* 314 */
        {cTan        , 488 }, /* 315 */
        {cPow        , 490 }, /* 316 */
        {cCosh       , 484 }, /* 317 */
        {cPow        , 493 }, /* 318 */
        {cTanh       , 488 }, /* 319 */
        {cAdd        , 500 }, /* 320 */
        {cTan        , 501 }, /* 321 */
        {cMul        , 502 }, /* 322 */
        {cTan        , 503 }, /* 323 */
        {cPow        , 508 }, /* 324 */
        {cAdd        , 509 }, /* 325 */
        {cAdd        , 510 }, /* 326 */
        {cPow        , 511 }, /* 327 */
        {cPow        , 514 }, /* 328 */
        {cMul        , 513 }, /* 329 */
        {cPow        , 516 }, /* 330 */
        {cAdd        , 518 }, /* 331 */
        {cAdd        , 519 }, /* 332 */
        {cPow        , 520 }, /* 333 */
        {cPow        , 522 }, /* 334 */
        {cSec        , 0 }, /* 335 */
        {cCsc        , 0 }, /* 336 */
        {cCot        , 0 }, /* 337 */
        {cRad        , 110 }, /* 338 */
        {cDeg        , 110 }, /* 339 */
    };

    const Rule rlist[] =
    {
        {1, ProduceNewTree,    6,	{ cAcosh      , 0 } }, /* 0 */
        {1, ProduceNewTree,    11,	{ cAsinh      , 0 } }, /* 1 */
        {1, ProduceNewTree,    19,	{ cAtanh      , 0 } }, /* 2 */
        {1, ProduceNewTree,    24,	{ cCosh       , 0 } }, /* 3 */
        {1, ProduceNewTree,    27,	{ cSinh       , 0 } }, /* 4 */
        {1, ProduceNewTree,    30,	{ cTan        , 0 } }, /* 5 */
        {1, ProduceNewTree,    35,	{ cTanh       , 0 } }, /* 6 */
        {1, ProduceNewTree,    37,	{ cNot        , 36 } }, /* 7 */
        {1, ProduceNewTree,    41,	{ cAbs        , 39 } }, /* 8 */
        {1, ProduceNewTree,    0,	{ cCos        , 42 } }, /* 9 */
        {1, ProduceNewTree,    47,	{ cCos        , 44 } }, /* 10 */
        {1, ProduceNewTree,    53,	{ cCos        , 50 } }, /* 11 */
        {1, ReplaceParams ,    52,	{ cCos        , 54 } }, /* 12 */
        {3, ProduceNewTree,    0,	{ cIf         , 55 } }, /* 13 */
        {3, ProduceNewTree,    60,	{ cIf         , 57 } }, /* 14 */
        {3, ProduceNewTree,    64,	{ cIf         , 61 } }, /* 15 */
        {3, ProduceNewTree,    68,	{ cIf         , 65 } }, /* 16 */
        {3, ProduceNewTree,    72,	{ cIf         , 69 } }, /* 17 */
        {3, ProduceNewTree,    78,	{ cIf         , 74 } }, /* 18 */
        {3, ProduceNewTree,    82,	{ cIf         , 79 } }, /* 19 */
        {3, ProduceNewTree,    86,	{ cIf         , 83 } }, /* 20 */
        {3, ProduceNewTree,    90,	{ cIf         , 87 } }, /* 21 */
        {3, ProduceNewTree,    94,	{ cIf         , 92 } }, /* 22 */
        {3, ProduceNewTree,    97,	{ cIf         , 95 } }, /* 23 */
        {1, ProduceNewTree,    102,	{ cLog        , 99 } }, /* 24 */
        {1, ProduceNewTree,    106,	{ cLog        , 104 } }, /* 25 */
        {1, ProduceNewTree,    114,	{ cLog        , 109 } }, /* 26 */
        {1, ProduceNewTree,    118,	{ cLog        , 116 } }, /* 27 */
        {2, ReplaceParams ,    0,	{ cMax        , 119 } }, /* 28 */
        {2, ReplaceParams ,    0,	{ cMin        , 119 } }, /* 29 */
        {2, ProduceNewTree,    0,	{ cPow        , 120 } }, /* 30 */
        {2, ReplaceParams ,    123,	{ cPow        , 122 } }, /* 31 */
        {2, ReplaceParams ,    125,	{ cPow        , 124 } }, /* 32 */
        {2, ReplaceParams ,    128,	{ cPow        , 127 } }, /* 33 */
        {2, ReplaceParams ,    123,	{ cPow        , 130 } }, /* 34 */
        {2, ReplaceParams ,    134,	{ cPow        , 132 } }, /* 35 */
        {2, ProduceNewTree,    0,	{ cPow        , 135 } }, /* 36 */
        {2, ReplaceParams ,    123,	{ cPow        , 139 } }, /* 37 */
        {2, ReplaceParams ,    141,	{ cPow        , 140 } }, /* 38 */
        {2, ReplaceParams ,    145,	{ cPow        , 143 } }, /* 39 */
        {2, ReplaceParams ,    145,	{ cPow        , 147 } }, /* 40 */
        {2, ReplaceParams ,    145,	{ cPow        , 148 } }, /* 41 */
        {2, ProduceNewTree,    153,	{ cPow        , 150 } }, /* 42 */
        {2, ReplaceParams ,    156,	{ cPow        , 155 } }, /* 43 */
        {1, ProduceNewTree,    0,	{ cSin        , 157 } }, /* 44 */
        {1, ProduceNewTree,    159,	{ cSin        , 44 } }, /* 45 */
        {1, ProduceNewTree,    160,	{ cSin        , 50 } }, /* 46 */
        {1, ProduceNewTree,    162,	{ cSin        , 54 } }, /* 47 */
        {2, ReplaceParams ,    167,	{ cAdd        , 165 } }, /* 48 */
        {2, ReplaceParams ,    170,	{ cAdd        , 169 } }, /* 49 */
        {2, ReplaceParams ,    174,	{ cAdd        , 172 } }, /* 50 */
        {2, ReplaceParams ,    176,	{ cAdd        , 175 } }, /* 51 */
        {2, ReplaceParams ,    182,	{ cAdd        , 178 } }, /* 52 */
        {2, ReplaceParams ,    195,	{ cAdd        , 185 } }, /* 53 */
        {2, ReplaceParams ,    208,	{ cAdd        , 198 } }, /* 54 */
        {2, ReplaceParams ,    213,	{ cAdd        , 211 } }, /* 55 */
        {2, ReplaceParams ,    219,	{ cAdd        , 215 } }, /* 56 */
        {2, ReplaceParams ,    223,	{ cAdd        , 222 } }, /* 57 */
        {2, ReplaceParams ,    226,	{ cAdd        , 225 } }, /* 58 */
        {2, ReplaceParams ,    232,	{ cAdd        , 228 } }, /* 59 */
        {2, ReplaceParams ,    236,	{ cAdd        , 234 } }, /* 60 */
        {2, ReplaceParams ,    239,	{ cAdd        , 237 } }, /* 61 */
        {2, ReplaceParams ,    243,	{ cAdd        , 241 } }, /* 62 */
        {2, ReplaceParams ,    232,	{ cAdd        , 245 } }, /* 63 */
        {2, ReplaceParams ,    249,	{ cAdd        , 246 } }, /* 64 */
        {2, ReplaceParams ,    252,	{ cAdd        , 250 } }, /* 65 */
        {3, ReplaceParams ,    261,	{ cAdd        , 256 } }, /* 66 */
        {2, ProduceNewTree,    266,	{ cMul        , 264 } }, /* 67 */
        {2, ReplaceParams ,    272,	{ cMul        , 269 } }, /* 68 */
        {2, ReplaceParams ,    276,	{ cMul        , 274 } }, /* 69 */
        {2, ReplaceParams ,    279,	{ cMul        , 277 } }, /* 70 */
        {2, ReplaceParams ,    286,	{ cMul        , 282 } }, /* 71 */
        {2, ReplaceParams ,    290,	{ cMul        , 287 } }, /* 72 */
        {2, ReplaceParams ,    294,	{ cMul        , 291 } }, /* 73 */
        {2, ReplaceParams ,    295,	{ cMul        , 250 } }, /* 74 */
        {2, ProduceNewTree,    176,	{ cEqual      , 296 } }, /* 75 */
        {2, ProduceNewTree,    297,	{ cNEqual     , 296 } }, /* 76 */
        {2, ProduceNewTree,    297,	{ cLess       , 296 } }, /* 77 */
        {2, ProduceNewTree,    176,	{ cLessOrEq   , 296 } }, /* 78 */
        {2, ProduceNewTree,    297,	{ cGreater    , 296 } }, /* 79 */
        {2, ProduceNewTree,    176,	{ cGreaterOrEq, 296 } }, /* 80 */
        {1, ProduceNewTree,    299,	{ cNot        , 298 } }, /* 81 */
        {1, ProduceNewTree,    298,	{ cNot        , 299 } }, /* 82 */
        {1, ProduceNewTree,    301,	{ cNot        , 300 } }, /* 83 */
        {1, ProduceNewTree,    303,	{ cNot        , 302 } }, /* 84 */
        {1, ProduceNewTree,    302,	{ cNot        , 303 } }, /* 85 */
        {1, ProduceNewTree,    300,	{ cNot        , 301 } }, /* 86 */
        {1, ProduceNewTree,    307,	{ cNot        , 305 } }, /* 87 */
        {1, ProduceNewTree,    309,	{ cNot        , 308 } }, /* 88 */
        {0, ReplaceParams ,    311,	{ cAnd        , 304 } }, /* 89 */
        {1, ReplaceParams ,    313,	{ cAnd        , 312 } }, /* 90 */
        {1, ReplaceParams ,    0,	{ cAnd        , 314 } }, /* 91 */
        {1, ReplaceParams ,    299,	{ cAnd        , 315 } }, /* 92 */
        {1, ReplaceParams ,    298,	{ cAnd        , 316 } }, /* 93 */
        {1, ReplaceParams ,    301,	{ cAnd        , 317 } }, /* 94 */
        {1, ReplaceParams ,    303,	{ cAnd        , 318 } }, /* 95 */
        {1, ReplaceParams ,    302,	{ cAnd        , 319 } }, /* 96 */
        {1, ReplaceParams ,    300,	{ cAnd        , 320 } }, /* 97 */
        {1, ReplaceParams ,    0,	{ cAnd        , 321 } }, /* 98 */
        {1, ReplaceParams ,    313,	{ cAnd        , 322 } }, /* 99 */
        {2, ReplaceParams ,    0,	{ cAnd        , 119 } }, /* 100 */
        {2, ProduceNewTree,    297,	{ cAnd        , 323 } }, /* 101 */
        {2, ProduceNewTree,    297,	{ cAnd        , 324 } }, /* 102 */
        {2, ProduceNewTree,    297,	{ cAnd        , 325 } }, /* 103 */
        {2, ProduceNewTree,    297,	{ cAnd        , 326 } }, /* 104 */
        {2, ProduceNewTree,    300,	{ cAnd        , 327 } }, /* 105 */
        {2, ProduceNewTree,    297,	{ cAnd        , 328 } }, /* 106 */
        {2, ProduceNewTree,    298,	{ cAnd        , 329 } }, /* 107 */
        {2, ProduceNewTree,    303,	{ cAnd        , 330 } }, /* 108 */
        {2, ReplaceParams ,    313,	{ cAnd        , 331 } }, /* 109 */
        {3, ReplaceParams ,    335,	{ cAnd        , 334 } }, /* 110 */
        {0, ReplaceParams ,    336,	{ cOr         , 304 } }, /* 111 */
        {1, ReplaceParams ,    313,	{ cOr         , 312 } }, /* 112 */
        {1, ReplaceParams ,    0,	{ cOr         , 314 } }, /* 113 */
        {1, ReplaceParams ,    299,	{ cOr         , 315 } }, /* 114 */
        {1, ReplaceParams ,    298,	{ cOr         , 316 } }, /* 115 */
        {1, ReplaceParams ,    301,	{ cOr         , 317 } }, /* 116 */
        {1, ReplaceParams ,    303,	{ cOr         , 318 } }, /* 117 */
        {1, ReplaceParams ,    302,	{ cOr         , 319 } }, /* 118 */
        {1, ReplaceParams ,    300,	{ cOr         , 320 } }, /* 119 */
        {1, ReplaceParams ,    0,	{ cOr         , 321 } }, /* 120 */
        {1, ReplaceParams ,    313,	{ cOr         , 322 } }, /* 121 */
        {2, ReplaceParams ,    0,	{ cOr         , 119 } }, /* 122 */
        {2, ProduceNewTree,    176,	{ cOr         , 323 } }, /* 123 */
        {2, ProduceNewTree,    176,	{ cOr         , 324 } }, /* 124 */
        {2, ProduceNewTree,    302,	{ cOr         , 337 } }, /* 125 */
        {2, ProduceNewTree,    299,	{ cOr         , 325 } }, /* 126 */
        {2, ProduceNewTree,    176,	{ cOr         , 326 } }, /* 127 */
        {2, ProduceNewTree,    302,	{ cOr         , 338 } }, /* 128 */
        {2, ProduceNewTree,    176,	{ cOr         , 328 } }, /* 129 */
        {2, ProduceNewTree,    301,	{ cOr         , 339 } }, /* 130 */
        {2, ProduceNewTree,    301,	{ cOr         , 340 } }, /* 131 */
        {2, ReplaceParams ,    313,	{ cOr         , 331 } }, /* 132 */
        {1, ProduceNewTree,    298,	{ cNotNot     , 298 } }, /* 133 */
        {1, ProduceNewTree,    299,	{ cNotNot     , 299 } }, /* 134 */
        {1, ProduceNewTree,    300,	{ cNotNot     , 300 } }, /* 135 */
        {1, ProduceNewTree,    302,	{ cNotNot     , 302 } }, /* 136 */
        {1, ProduceNewTree,    303,	{ cNotNot     , 303 } }, /* 137 */
        {1, ProduceNewTree,    301,	{ cNotNot     , 301 } }, /* 138 */
        {1, ProduceNewTree,    341,	{ cNotNot     , 36 } }, /* 139 */
        {1, ProduceNewTree,    345,	{ cNotNot     , 343 } }, /* 140 */
        {1, ProduceNewTree,    347,	{ cNotNot     , 346 } }, /* 141 */
        {1, ReplaceParams ,    0,	{ cNotNot     , 37 } }, /* 142 */
        {2, ReplaceParams ,    349,	{ cAtan2      , 348 } }, /* 143 */
        {1, ProduceNewTree,    351,	{ cCosh       , 350 } }, /* 144 */
        {1, ReplaceParams ,    110,	{ cCosh       , 353 } }, /* 145 */
        {3, ProduceNewTree,    0,	{ cIf         , 55 } }, /* 146 */
        {3, ProduceNewTree,    60,	{ cIf         , 57 } }, /* 147 */
        {3, ProduceNewTree,    64,	{ cIf         , 61 } }, /* 148 */
        {3, ProduceNewTree,    68,	{ cIf         , 65 } }, /* 149 */
        {3, ProduceNewTree,    72,	{ cIf         , 69 } }, /* 150 */
        {3, ProduceNewTree,    78,	{ cIf         , 74 } }, /* 151 */
        {3, ProduceNewTree,    82,	{ cIf         , 79 } }, /* 152 */
        {3, ProduceNewTree,    86,	{ cIf         , 83 } }, /* 153 */
        {3, ProduceNewTree,    90,	{ cIf         , 87 } }, /* 154 */
        {3, ProduceNewTree,    94,	{ cIf         , 92 } }, /* 155 */
        {3, ProduceNewTree,    97,	{ cIf         , 95 } }, /* 156 */
        {1, ProduceNewTree,    102,	{ cLog        , 99 } }, /* 157 */
        {1, ProduceNewTree,    106,	{ cLog        , 104 } }, /* 158 */
        {1, ProduceNewTree,    358,	{ cLog        , 357 } }, /* 159 */
        {1, ProduceNewTree,    350,	{ cLog        , 362 } }, /* 160 */
        {1, ProduceNewTree,    114,	{ cLog        , 109 } }, /* 161 */
        {1, ProduceNewTree,    118,	{ cLog        , 116 } }, /* 162 */
        {1, ProduceNewTree,    370,	{ cLog        , 368 } }, /* 163 */
        {1, ProduceNewTree,    380,	{ cLog        , 376 } }, /* 164 */
        {2, ReplaceParams ,    0,	{ cMax        , 119 } }, /* 165 */
        {2, ReplaceParams ,    0,	{ cMin        , 119 } }, /* 166 */
        {2, ProduceNewTree,    0,	{ cPow        , 120 } }, /* 167 */
        {2, ReplaceParams ,    123,	{ cPow        , 122 } }, /* 168 */
        {2, ReplaceParams ,    125,	{ cPow        , 124 } }, /* 169 */
        {2, ReplaceParams ,    128,	{ cPow        , 127 } }, /* 170 */
        {2, ReplaceParams ,    123,	{ cPow        , 130 } }, /* 171 */
        {2, ReplaceParams ,    134,	{ cPow        , 132 } }, /* 172 */
        {2, ProduceNewTree,    0,	{ cPow        , 135 } }, /* 173 */
        {2, ReplaceParams ,    123,	{ cPow        , 139 } }, /* 174 */
        {2, ReplaceParams ,    141,	{ cPow        , 140 } }, /* 175 */
        {2, ReplaceParams ,    145,	{ cPow        , 143 } }, /* 176 */
        {2, ReplaceParams ,    145,	{ cPow        , 147 } }, /* 177 */
        {2, ReplaceParams ,    145,	{ cPow        , 148 } }, /* 178 */
        {2, ProduceNewTree,    153,	{ cPow        , 150 } }, /* 179 */
        {2, ReplaceParams ,    156,	{ cPow        , 155 } }, /* 180 */
        {1, ProduceNewTree,    381,	{ cSinh       , 358 } }, /* 181 */
        {1, ProduceNewTree,    383,	{ cSinh       , 353 } }, /* 182 */
        {1, ProduceNewTree,    0,	{ cTan        , 384 } }, /* 183 */
        {1, ProduceNewTree,    386,	{ cTan        , 353 } }, /* 184 */
        {1, ProduceNewTree,    388,	{ cTanh       , 353 } }, /* 185 */
        {2, ProduceNewTree,    390,	{ cAdd        , 356 } }, /* 186 */
        {2, ProduceNewTree,    392,	{ cAdd        , 361 } }, /* 187 */
        {2, ReplaceParams ,    398,	{ cAdd        , 393 } }, /* 188 */
        {2, ReplaceParams ,    401,	{ cAdd        , 399 } }, /* 189 */
        {2, ReplaceParams ,    404,	{ cAdd        , 402 } }, /* 190 */
        {2, ReplaceParams ,    407,	{ cAdd        , 406 } }, /* 191 */
        {2, ReplaceParams ,    410,	{ cAdd        , 409 } }, /* 192 */
        {2, ReplaceParams ,    174,	{ cAdd        , 172 } }, /* 193 */
        {2, ReplaceParams ,    417,	{ cAdd        , 413 } }, /* 194 */
        {2, ReplaceParams ,    421,	{ cAdd        , 419 } }, /* 195 */
        {2, ReplaceParams ,    423,	{ cAdd        , 422 } }, /* 196 */
        {2, ReplaceParams ,    427,	{ cAdd        , 425 } }, /* 197 */
        {2, ReplaceParams ,    430,	{ cAdd        , 428 } }, /* 198 */
        {2, ReplaceParams ,    433,	{ cAdd        , 431 } }, /* 199 */
        {2, ReplaceParams ,    437,	{ cAdd        , 435 } }, /* 200 */
        {2, ReplaceParams ,    440,	{ cAdd        , 438 } }, /* 201 */
        {2, ReplaceParams ,    445,	{ cAdd        , 443 } }, /* 202 */
        {2, ReplaceParams ,    449,	{ cAdd        , 447 } }, /* 203 */
        {2, ReplaceParams ,    249,	{ cAdd        , 246 } }, /* 204 */
        {2, ReplaceParams ,    252,	{ cAdd        , 250 } }, /* 205 */
        {2, ProduceNewTree,    266,	{ cMul        , 264 } }, /* 206 */
        {2, ReplaceParams ,    272,	{ cMul        , 269 } }, /* 207 */
        {2, ReplaceParams ,    53,	{ cMul        , 451 } }, /* 208 */
        {2, ReplaceParams ,    453,	{ cMul        , 452 } }, /* 209 */
        {2, ReplaceParams ,    455,	{ cMul        , 454 } }, /* 210 */
        {2, ReplaceParams ,    457,	{ cMul        , 456 } }, /* 211 */
        {2, ReplaceParams ,    276,	{ cMul        , 274 } }, /* 212 */
        {2, ReplaceParams ,    279,	{ cMul        , 277 } }, /* 213 */
        {2, ReplaceParams ,    461,	{ cMul        , 459 } }, /* 214 */
        {2, ReplaceParams ,    465,	{ cMul        , 463 } }, /* 215 */
        {2, ReplaceParams ,    467,	{ cMul        , 466 } }, /* 216 */
        {2, ReplaceParams ,    423,	{ cMul        , 468 } }, /* 217 */
        {2, ReplaceParams ,    480,	{ cMul        , 474 } }, /* 218 */
        {2, ReplaceParams ,    482,	{ cMul        , 481 } }, /* 219 */
        {2, ReplaceParams ,    489,	{ cMul        , 486 } }, /* 220 */
        {2, ReplaceParams ,    492,	{ cMul        , 491 } }, /* 221 */
        {2, ReplaceParams ,    495,	{ cMul        , 494 } }, /* 222 */
        {2, ReplaceParams ,    497,	{ cMul        , 496 } }, /* 223 */
        {2, ReplaceParams ,    499,	{ cMul        , 498 } }, /* 224 */
        {2, ReplaceParams ,    176,	{ cMul        , 504 } }, /* 225 */
        {2, ReplaceParams ,    407,	{ cMul        , 505 } }, /* 226 */
        {2, ReplaceParams ,    507,	{ cMul        , 506 } }, /* 227 */
        {2, ReplaceParams ,    513,	{ cMul        , 512 } }, /* 228 */
        {2, ReplaceParams ,    517,	{ cMul        , 515 } }, /* 229 */
        {2, ReplaceParams ,    492,	{ cMul        , 521 } }, /* 230 */
        {2, ReplaceParams ,    465,	{ cMul        , 523 } }, /* 231 */
        {2, ReplaceParams ,    294,	{ cMul        , 291 } }, /* 232 */
        {2, ReplaceParams ,    295,	{ cMul        , 250 } }, /* 233 */
        {2, ProduceNewTree,    176,	{ cEqual      , 296 } }, /* 234 */
        {2, ProduceNewTree,    297,	{ cNEqual     , 296 } }, /* 235 */
        {2, ProduceNewTree,    297,	{ cLess       , 296 } }, /* 236 */
        {2, ProduceNewTree,    176,	{ cLessOrEq   , 296 } }, /* 237 */
        {2, ProduceNewTree,    297,	{ cGreater    , 296 } }, /* 238 */
        {2, ProduceNewTree,    176,	{ cGreaterOrEq, 296 } }, /* 239 */
        {1, ProduceNewTree,    299,	{ cNot        , 298 } }, /* 240 */
        {1, ProduceNewTree,    298,	{ cNot        , 299 } }, /* 241 */
        {1, ProduceNewTree,    301,	{ cNot        , 300 } }, /* 242 */
        {1, ProduceNewTree,    303,	{ cNot        , 302 } }, /* 243 */
        {1, ProduceNewTree,    302,	{ cNot        , 303 } }, /* 244 */
        {1, ProduceNewTree,    300,	{ cNot        , 301 } }, /* 245 */
        {1, ProduceNewTree,    307,	{ cNot        , 305 } }, /* 246 */
        {1, ProduceNewTree,    309,	{ cNot        , 308 } }, /* 247 */
        {0, ReplaceParams ,    311,	{ cAnd        , 304 } }, /* 248 */
        {1, ReplaceParams ,    313,	{ cAnd        , 312 } }, /* 249 */
        {1, ReplaceParams ,    0,	{ cAnd        , 314 } }, /* 250 */
        {1, ReplaceParams ,    299,	{ cAnd        , 315 } }, /* 251 */
        {1, ReplaceParams ,    298,	{ cAnd        , 316 } }, /* 252 */
        {1, ReplaceParams ,    301,	{ cAnd        , 317 } }, /* 253 */
        {1, ReplaceParams ,    303,	{ cAnd        , 318 } }, /* 254 */
        {1, ReplaceParams ,    302,	{ cAnd        , 319 } }, /* 255 */
        {1, ReplaceParams ,    300,	{ cAnd        , 320 } }, /* 256 */
        {1, ReplaceParams ,    0,	{ cAnd        , 321 } }, /* 257 */
        {1, ReplaceParams ,    313,	{ cAnd        , 322 } }, /* 258 */
        {2, ReplaceParams ,    0,	{ cAnd        , 119 } }, /* 259 */
        {2, ProduceNewTree,    297,	{ cAnd        , 323 } }, /* 260 */
        {2, ProduceNewTree,    297,	{ cAnd        , 324 } }, /* 261 */
        {2, ProduceNewTree,    297,	{ cAnd        , 325 } }, /* 262 */
        {2, ProduceNewTree,    297,	{ cAnd        , 326 } }, /* 263 */
        {2, ProduceNewTree,    300,	{ cAnd        , 327 } }, /* 264 */
        {2, ProduceNewTree,    297,	{ cAnd        , 328 } }, /* 265 */
        {2, ProduceNewTree,    298,	{ cAnd        , 329 } }, /* 266 */
        {2, ProduceNewTree,    303,	{ cAnd        , 330 } }, /* 267 */
        {2, ReplaceParams ,    313,	{ cAnd        , 331 } }, /* 268 */
        {3, ReplaceParams ,    335,	{ cAnd        , 334 } }, /* 269 */
        {0, ReplaceParams ,    336,	{ cOr         , 304 } }, /* 270 */
        {1, ReplaceParams ,    313,	{ cOr         , 312 } }, /* 271 */
        {1, ReplaceParams ,    0,	{ cOr         , 314 } }, /* 272 */
        {1, ReplaceParams ,    299,	{ cOr         , 315 } }, /* 273 */
        {1, ReplaceParams ,    298,	{ cOr         , 316 } }, /* 274 */
        {1, ReplaceParams ,    301,	{ cOr         , 317 } }, /* 275 */
        {1, ReplaceParams ,    303,	{ cOr         , 318 } }, /* 276 */
        {1, ReplaceParams ,    302,	{ cOr         , 319 } }, /* 277 */
        {1, ReplaceParams ,    300,	{ cOr         , 320 } }, /* 278 */
        {1, ReplaceParams ,    0,	{ cOr         , 321 } }, /* 279 */
        {1, ReplaceParams ,    313,	{ cOr         , 322 } }, /* 280 */
        {2, ReplaceParams ,    0,	{ cOr         , 119 } }, /* 281 */
        {2, ProduceNewTree,    176,	{ cOr         , 323 } }, /* 282 */
        {2, ProduceNewTree,    176,	{ cOr         , 324 } }, /* 283 */
        {2, ProduceNewTree,    302,	{ cOr         , 337 } }, /* 284 */
        {2, ProduceNewTree,    299,	{ cOr         , 325 } }, /* 285 */
        {2, ProduceNewTree,    176,	{ cOr         , 326 } }, /* 286 */
        {2, ProduceNewTree,    302,	{ cOr         , 338 } }, /* 287 */
        {2, ProduceNewTree,    176,	{ cOr         , 328 } }, /* 288 */
        {2, ProduceNewTree,    301,	{ cOr         , 339 } }, /* 289 */
        {2, ProduceNewTree,    301,	{ cOr         , 340 } }, /* 290 */
        {2, ReplaceParams ,    313,	{ cOr         , 331 } }, /* 291 */
        {1, ProduceNewTree,    298,	{ cNotNot     , 298 } }, /* 292 */
        {1, ProduceNewTree,    299,	{ cNotNot     , 299 } }, /* 293 */
        {1, ProduceNewTree,    300,	{ cNotNot     , 300 } }, /* 294 */
        {1, ProduceNewTree,    302,	{ cNotNot     , 302 } }, /* 295 */
        {1, ProduceNewTree,    303,	{ cNotNot     , 303 } }, /* 296 */
        {1, ProduceNewTree,    301,	{ cNotNot     , 301 } }, /* 297 */
        {1, ProduceNewTree,    341,	{ cNotNot     , 36 } }, /* 298 */
        {1, ProduceNewTree,    345,	{ cNotNot     , 343 } }, /* 299 */
        {1, ProduceNewTree,    347,	{ cNotNot     , 346 } }, /* 300 */
        {1, ReplaceParams ,    0,	{ cNotNot     , 37 } }, /* 301 */
        {2, ReplaceParams ,    525,	{ cPow        , 524 } }, /* 302 */
        {2, ReplaceParams ,    527,	{ cPow        , 526 } }, /* 303 */
        {2, ReplaceParams ,    529,	{ cPow        , 528 } }, /* 304 */
        {1, ProduceNewTree,    531,	{ cMul        , 530 } }, /* 305 */
        {1, ProduceNewTree,    533,	{ cMul        , 532 } }, /* 306 */
    };
}

namespace FPoptimizer_Grammar
{
    const GrammarPack pack =
    {
        clist, plist, mlist, flist, rlist,
        {
            {0, 8 }, /* 0 */
            {8, 135 }, /* 1 */
            {143, 159 }, /* 2 */
            {302, 5 }, /* 3 */
        }
    };
}

#include <algorithm>
#include <cmath>
#include <map>
#include <assert.h>

#include "io/fparser/fpconfig.h"
#include "io/fparser/fparser.h"
#include "io/fparser/fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;

//#define DEBUG_SUBSTITUTIONS

#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
namespace FPoptimizer_Grammar
{
    void DumpTree(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o = std::cout);
    void DumpHashes(const FPoptimizer_CodeTree::CodeTree& tree);
}
#endif

namespace
{
    /* I have heard that std::equal_range() is practically worthless
     * due to the insane limitation that the two parameters for Comp() must
     * be of the same type. Hence we must reinvent the wheel and implement
     * our own here. This is practically identical to the one from
     * GNU libstdc++, except rewritten. -Bisqwit
     */
    template<typename It, typename T, typename Comp>
    std::pair<It, It>
    MyEqualRange(It first, It last, const T& val, Comp comp)
    {
        size_t len = last-first;
        while(len > 0)
        {
            size_t half = len/2;
            It middle(first); middle += half;
            if(comp(*middle, val))
            {
                first = middle;
                ++first;
                len = len - half - 1;
            }
            else if(comp(val, *middle))
            {
                len = half;
            }
            else
            {
                // The following implements this:
                // // left = lower_bound(first, middle, val, comp);
                It left(first);
              {///
                It& first2 = left;
                It last2(middle);
                size_t len2 = last2-first2;
                while(len2 > 0)
                {
                    size_t half2 = len2 / 2;
                    It middle2(first2); middle2 += half2;
                    if(comp(*middle2, val))
                    {
                        first2 = middle2;
                        ++first2;
                        len2 = len2 - half2 - 1;
                    }
                    else
                        len2 = half2;
                }
                // left = first2;  - not needed, already happens due to reference
              }///
                first += len;
                // The following implements this:
                // // right = upper_bound(++middle, first, val, comp);
                It right(++middle);
              {///
                It& first2 = right;
                It& last2 = first;
                size_t len2 = last2-first2;
                while(len2 > 0)
                {
                    size_t half2 = len2 / 2;
                    It middle2(first2); middle2 += half2;
                    if(comp(val, *middle2))
                        len2 = half2;
                    else
                    {
                        first2 = middle2;
                        ++first2;
                        len2 = len2 - half2 - 1;
                    }
                }
                // right = first2;  - not needed, already happens due to reference
              }///
                return std::pair<It,It> (left,right);
            }
        }
        return std::pair<It,It> (first,first);
    }
}

namespace FPoptimizer_CodeTree
{
    void CodeTree::ConstantFolding()
    {
        using namespace std;

        // Insert here any hardcoded constant-folding optimizations
        // that you want to be done whenever a new subtree is generated.
        /* Not recursive. */

        double const_value = 1.0;
        size_t which_param = 0;

        if(Opcode != cImmed)
        {
            MinMaxTree p = CalculateResultBoundaries();
            if(p.has_min && p.has_max && p.min == p.max)
            {
                // Replace us with this immed
                const_value = p.min;
                goto ReplaceTreeWithConstValue;
            }
        }


        /* Sub-list assimilation prepass */
        switch( (OPCODE) Opcode)
        {
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            {
                /* If the list contains another list of the same kind, assimilate it */
                for(size_t a=Params.size(); a-- > 0; )
                    if(Params[a].param->Opcode == Opcode)
                    {
                        // Assimilate its children and remove it
                        CodeTreeP tree = Params[a].param;
                        bool had_sign = Params[a].sign;

                        if(Opcode == cMul && had_sign)
                        {
                            int n_divisions_old = 1;
                            int n_divisions_new = 0;

                            for(size_t b=0; b<tree->Params.size(); ++b)
                                if(tree->Params[b].sign)
                                    ++n_divisions_old;
                                else
                                    ++n_divisions_new;

                            if(n_divisions_new > n_divisions_old)
                            {
                                // Don't assimilate if the number of divisions increases
                                // in the process.  x/(y*z) is better than x/y/z
                                continue;
                            }
                        }

                        Params.erase(Params.begin()+a);
                        for(size_t b=0; b<tree->Params.size(); ++b)
                            AddParam( Param(tree->Params[b].param,
                                            had_sign ^ tree->Params[b].sign) );
                    }
                break;
            }
            case cAnd:
            case cOr:
            {
                /* If the list contains another list of the same kind, assimilate it */
                for(size_t a=Params.size(); a-- > 0; )
                    if(Params[a].param->Opcode == Opcode
                    && Params[a].sign == false)
                    {
                        // Assimilate its children and remove it
                        CodeTreeP tree = Params[a].param;
                        Params.erase(Params.begin()+a);
                        for(size_t b=0; b<tree->Params.size(); ++b)
                            AddParam(tree->Params[b]);
                    }
                break;
            }
            default: break;
        }

        /* Constant folding */
        switch( (OPCODE) Opcode)
        {
            case cEval:
            case cImmed:
                break; // nothing to do
            case cVar:
                break; // nothing to do

            ReplaceTreeWithOne:
                const_value = 1.0;
                goto ReplaceTreeWithConstValue;
            ReplaceTreeWithZero:
                const_value = 0.0;
            ReplaceTreeWithConstValue:
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Replacing "; FPoptimizer_Grammar::DumpTree(*this);
                std::cout << " with const value " << const_value << "\n";
              #endif
                Params.clear();
                Opcode = cImmed;
                Value  = const_value;
                break;
            ReplaceTreeWithParam0:
                which_param = 0;
            ReplaceTreeWithParam:
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before replace: "; FPoptimizer_Grammar::DumpTree(*this);
                std::cout << "\n";
              #endif
                Opcode = Params[which_param].param->Opcode;
                Var    = Params[which_param].param->Var;
                Value  = Params[which_param].param->Value;
                Params.swap(Params[which_param].param->Params);
                for(size_t a=0; a<Params.size(); ++a)
                    Params[a].param->Parent = this;
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After replace: "; FPoptimizer_Grammar::DumpTree(*this);
                std::cout << "\n";
              #endif
                break;

            case cAnd:
            {
                // If the and-list contains an expression that evaluates to approx. zero,
                // the whole list evaluates to zero.
                // If all expressions within the and-list evaluate to approx. nonzero,
                // the whole list evaluates to one.
                bool all_values_are_nonzero = true;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    MinMaxTree p = Params[a].param->CalculateResultBoundaries();
                    if(p.has_min && p.has_max
                    && p.min > -0.5 && p.max < 0.5) // -0.5 < x < 0.5 = zero
                    {
                        if(!Params[a].sign) goto ReplaceTreeWithZero;
                        all_values_are_nonzero = false;
                    }
                    else if( (p.has_max && p.max <= -0.5)
                          || (p.has_min && p.min >= 0.5)) // |x| >= 0.5  = nonzero
                    {
                        if(Params[a].sign) goto ReplaceTreeWithZero;
                    }
                    else
                        all_values_are_nonzero = false;
                }
                if(all_values_are_nonzero) goto ReplaceTreeWithOne;
                if(Params.size() == 1 && !Params[0].sign)
                {
                    // Replace self with the single operand
                    Opcode = Params[0].sign ? cNot : cNotNot;
                    Params[0].sign = false;
                }
                if(Params.empty()) goto ReplaceTreeWithZero;
                break;
            }
            case cOr:
            {
                // If the or-list contains an expression that evaluates to approx. nonzero,
                // the whole list evaluates to one.
                // If all expressions within the and-list evaluate to approx. zero,
                // the whole list evaluates to zero.
                bool all_values_are_zero = true;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    MinMaxTree p = Params[a].param->CalculateResultBoundaries();
                    if(p.has_min && p.has_max
                    && p.min > -0.5 && p.max < 0.5) // -0.5 < x < 0.5 = zero
                    {
                        if(Params[a].sign) goto ReplaceTreeWithOne;
                    }
                    else if( (p.has_max && p.max <= -0.5)
                          || (p.has_min && p.min >= 0.5)) // |x| >= 0.5  = nonzero
                    {
                        if(!Params[a].sign) goto ReplaceTreeWithOne;
                        all_values_are_zero = false;
                    }
                    else
                        all_values_are_zero = false;
                }
                if(all_values_are_zero) goto ReplaceTreeWithZero;
                if(Params.size() == 1)
                {
                    // Replace self with the single operand
                    Opcode = Params[0].sign ? cNot : cNotNot;
                    Params[0].sign = false;
                }
                if(Params.empty()) goto ReplaceTreeWithOne;
                break;
            }
            case cNot:
            {
                // If the sub-expression evaluates to approx. zero, yield one.
                // If the sub-expression evaluates to approx. nonzero, yield zero.
                MinMaxTree p = Params[0].param->CalculateResultBoundaries();
                if(p.has_min && p.has_max
                && p.min > -0.5 && p.max < 0.5) // -0.5 < x < 0.5 = zero
                {
                    goto ReplaceTreeWithOne;
                }
                else if( (p.has_max && p.max <= -0.5)
                      || (p.has_min && p.min >= 0.5)) // |x| >= 0.5  = nonzero
                    goto ReplaceTreeWithZero;
                break;
            }
            case cNotNot:
            {
                // If the sub-expression evaluates to approx. zero, yield zero.
                // If the sub-expression evaluates to approx. nonzero, yield one.
                MinMaxTree p = Params[0].param->CalculateResultBoundaries();
                if(p.has_min && p.has_max
                && p.min > -0.5 && p.max < 0.5) // -0.5 < x < 0.5 = zero
                {
                    goto ReplaceTreeWithZero;
                }
                else if( (p.has_max && p.max <= -0.5)
                      || (p.has_min && p.min >= 0.5)) // |x| >= 0.5  = nonzero
                    goto ReplaceTreeWithOne;
                break;
            }
            case cIf:
            {
                // If the sub-expression evaluates to approx. zero, yield param3.
                // If the sub-expression evaluates to approx. nonzero, yield param2.
                MinMaxTree p = Params[0].param->CalculateResultBoundaries();
                if(p.has_min && p.has_max
                && p.min > -0.5 && p.max < 0.5) // -0.5 < x < 0.5 = zero
                {
                    which_param = 2;
                    goto ReplaceTreeWithParam;
                }
                else if( (p.has_max && p.max <= -0.5)
                      || (p.has_min && p.min >= 0.5)) // |x| >= 0.5  = nonzero
                {
                    which_param = 1;
                    goto ReplaceTreeWithParam;
                }
                break;
            }
            case cMul:
            {
            NowWeAreMulGroup: ;
                // If one sub-expression evalutes to exact zero, yield zero.
                double mul_immed_sum = 1.0;
                size_t n_mul_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    if(!Params[a].param->IsImmed()) continue;
                    // ^ Only check constant values
                    double immed = Params[a].param->GetImmed();
                    if(FloatEqual(immed, 0.0)) goto ReplaceTreeWithZero;
                    if(FloatEqual(immed, 1.0)) needs_resynth = true;
                    if(Params[a].sign) immed = 1.0 / immed;
                    mul_immed_sum *= immed; ++n_mul_immeds;
                }
                // Merge immeds.
                if(n_mul_immeds > 1) needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                    //std::cout << "cMul: Will add new immed " << mul_immed_sum << "\n";
                    for(size_t a=Params.size(); a-->0; )
                        if(Params[a].param->IsImmed())
                        {
                            //std::cout << " - For that, deleting immed " << Params[a].param->GetImmed();
                            //std::cout << "\n";
                            Params.erase(Params.begin()+a);
                        }
                    if(!FloatEqual(mul_immed_sum, 1.0))
                        AddParam( Param(new CodeTree(mul_immed_sum), false) );
                }
                if(Params.size() == 1 && !Params[0].sign)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                if(Params.empty()) goto ReplaceTreeWithOne;
                break;
            }
            case cAdd:
            {
                double immed_sum = 0.0;
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    if(!Params[a].param->IsImmed()) continue;
                    // ^ Only check constant values
                    double immed = Params[a].param->GetImmed();
                    if(FloatEqual(immed, 0.0)) needs_resynth = true;
                    if(Params[a].sign) immed = -immed;
                    immed_sum += immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1) needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                    //std::cout << "cAdd: Will add new immed " << immed_sum << "\n";
                    //std::cout << "In: "; FPoptimizer_Grammar::DumpTree(*this);
                    //std::cout << "\n";

                    for(size_t a=Params.size(); a-->0; )
                        if(Params[a].param->IsImmed())
                        {
                            //std::cout << " - For that, deleting immed " << Params[a].param->GetImmed();
                            //std::cout << "\n";
                            Params.erase(Params.begin()+a);
                        }
                    if(!FloatEqual(immed_sum, 0.0))
                        AddParam( Param(new CodeTree(immed_sum), false) );
                }
                if(Params.size() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                if(Params.empty()) goto ReplaceTreeWithZero;
                break;
            }
            case cMin:
            {
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the larger one.
                 *
                 * Algorithm: 1. figure out the smallest maximum of all operands.
                 *            2. eliminate all operands where their minimum is
                 *               larger than the selected maximum.
                 */
                MinMaxTree smallest_maximum;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    MinMaxTree p = Params[a].param->CalculateResultBoundaries();
                    if(p.has_max && (!smallest_maximum.has_max || p.max < smallest_maximum.max))
                    {
                        smallest_maximum.max = p.max;
                        smallest_maximum.has_max = true;
                }   }
                if(smallest_maximum.has_max)
                    for(size_t a=Params.size(); a-- > 0; )
                    {
                        MinMaxTree p = Params[a].param->CalculateResultBoundaries();
                        if(p.has_min && p.min > smallest_maximum.max)
                            Params.erase(Params.begin() + a);
                    }
                //fprintf(stderr, "Remains: %u\n", (unsigned)Params.size());
                if(Params.size() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }
            case cMax:
            {
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the smaller one.
                 *
                 * Algorithm: 1. figure out the biggest minimum of all operands.
                 *            2. eliminate all operands where their maximum is
                 *               smaller than the selected minimum.
                 */
                MinMaxTree biggest_minimum;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    MinMaxTree p = Params[a].param->CalculateResultBoundaries();
                    if(p.has_min && (!biggest_minimum.has_min || p.min > biggest_minimum.min))
                    {
                        biggest_minimum.min = p.min;
                        biggest_minimum.has_min = true;
                }   }
                if(biggest_minimum.has_min)
                {
                    //fprintf(stderr, "Removing all where max < %g\n", biggest_minimum.min);
                    for(size_t a=Params.size(); a-- > 0; )
                    {
                        MinMaxTree p = Params[a].param->CalculateResultBoundaries();
                        if(p.has_max && p.max < biggest_minimum.min)
                        {
                            //fprintf(stderr, "Removing %g\n", p.max);
                            Params.erase(Params.begin() + a);
                        }
                    }
                }
                //fprintf(stderr, "Remains: %u\n", (unsigned)Params.size());
                if(Params.size() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }

            case cEqual:
            {
                /* If we know the two operands' ranges don't overlap, we get zero.
                 * The opposite is more complex and is done in .dat code.
                 */
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if((p0.has_max && p1.has_min && p1.min > p0.max)
                || (p1.has_max && p0.has_min && p0.min > p1.max))
                    goto ReplaceTreeWithZero;
                break;
            }

            case cNEqual:
            {
                /* If we know the two operands' ranges don't overlap, we get one.
                 * The opposite is more complex and is done in .dat code.
                 */
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if((p0.has_max && p1.has_min && p1.min > p0.max)
                || (p1.has_max && p0.has_min && p0.min > p1.max))
                    goto ReplaceTreeWithOne;
                break;
            }

            case cLess:
            {
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if(p0.has_max && p1.has_min && p0.max < p1.min)
                    goto ReplaceTreeWithOne; // We know p0 < p1
                if(p1.has_max && p0.has_min && p1.max <= p0.min)
                    goto ReplaceTreeWithZero; // We know p1 >= p0
                break;
            }

            case cLessOrEq:
            {
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if(p0.has_max && p1.has_min && p0.max <= p1.min)
                    goto ReplaceTreeWithOne; // We know p0 <= p1
                if(p1.has_max && p0.has_min && p1.max < p0.min)
                    goto ReplaceTreeWithZero; // We know p1 > p0
                break;
            }

            case cGreater:
            {
                // Note: Eq case not handled
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if(p0.has_max && p1.has_min && p0.max <= p1.min)
                    goto ReplaceTreeWithZero; // We know p0 <= p1
                if(p1.has_max && p0.has_min && p1.max < p0.min)
                    goto ReplaceTreeWithOne; // We know p1 > p0
                break;
            }

            case cGreaterOrEq:
            {
                // Note: Eq case not handled
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if(p0.has_max && p1.has_min && p0.max < p1.min)
                    goto ReplaceTreeWithZero; // We know p0 < p1
                if(p1.has_max && p0.has_min && p1.max <= p0.min)
                    goto ReplaceTreeWithOne; // We know p1 >= p0
                break;
            }

            case cAbs:
            {
                /* If we know the operand is always positive, cAbs is redundant.
                 * If we know the operand is always negative, use actual negation.
                 */
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                if(p0.has_min && p0.min >= 0.0)
                    goto ReplaceTreeWithParam0;
                if(p0.has_max && p0.max <= NEGATIVE_MAXIMUM)
                {
                    /* abs(negative) = negative*-1 */
                    Opcode = cMul;
                    AddParam( Param(new CodeTree(-1.0), false) );
                    /* The caller of ConstantFolding() will do Sort() and Rehash() next.
                     * Thus, no need to do it here. */
                    /* We were changed into a cMul group. Do cMul folding. */
                    goto NowWeAreMulGroup;
                }
                /* If the operand is a cMul group, find elements
                 * that are always positive and always negative,
                 * and move them out, e.g. abs(p*n*x*y) = p*(-n)*abs(x*y)
                 */
                if(Params[0].param->Opcode == cMul)
                {
                    CodeTree& p = *Params[0].param;
                    std::vector<Param> pos_set;
                    std::vector<Param> neg_set;
                    for(size_t a=0; a<p.Params.size(); ++a)
                    {
                        p0 = p.Params[a].param->CalculateResultBoundaries();
                        if(p0.has_min && p0.min >= 0.0)
                            { pos_set.push_back(p.Params[a]); }
                        if(p0.has_max && p0.max <= NEGATIVE_MAXIMUM)
                            { neg_set.push_back(p.Params[a]); }
                    }
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "Abs: mul group has " << pos_set.size()
                              << " pos, " << neg_set.size() << "neg\n";
                #endif
                    if(!pos_set.empty() || !neg_set.empty())
                    {
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-Before: ";
                        FPoptimizer_Grammar::DumpTree(*this);
                        std::cout << "\n" << std::flush;
                        FPoptimizer_Grammar::DumpHashes(*this);
                #endif
                        for(size_t a=p.Params.size(); a-- > 0; )
                        {
                            p0 = p.Params[a].param->CalculateResultBoundaries();
                            if((p0.has_min && p0.min >= 0.0)
                            || (p0.has_max && p0.max <= NEGATIVE_MAXIMUM))
                                p.Params.erase(p.Params.begin() + a);

                            /* Here, p*n*x*y -> x*y.
                             * p is saved in pos_set[]
                             * n is saved in neg_set[]
                             */
                        }
                        p.ConstantFolding();
                        p.Sort();

                        CodeTreeP subtree = new CodeTree;
                        p.Parent = &*subtree;
                        subtree->Opcode = cAbs;
                        subtree->Params.swap(Params);
                        subtree->ConstantFolding();
                        subtree->Sort();
                        subtree->Rehash(false); // hash it and its children.

                        /* Now:
                         * subtree = Abs(x*y)
                         * this    = Abs()
                         */

                        Opcode = cMul;
                        for(size_t a=0; a<pos_set.size(); ++a)
                            AddParam(pos_set[a]);
                        AddParam(Param(subtree, false));
                        /* Now:
                         * this    = p * Abs(x*y)
                         */
                        if(!neg_set.empty())
                        {
                            for(size_t a=0; a<neg_set.size(); ++a)
                                AddParam(neg_set[a]);
                            AddParam( Param(new CodeTree(-1.0), false) );
                            /* Now:
                             * this = p * n * -1 * Abs(x*y)
                             */
                        }
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-After: ";
                        FPoptimizer_Grammar::DumpTree(*this);
                        std::cout << "\n" << std::flush;
                        FPoptimizer_Grammar::DumpHashes(*this);
                #endif
                        /* We were changed into a cMul group. Do cMul folding. */
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            }

            #define HANDLE_UNARY_CONST_FUNC(funcname) \
                if(Params[0].param->IsImmed()) \
                    { const_value = funcname(Params[0].param->GetImmed()); \
                      goto ReplaceTreeWithConstValue; }

            case cLog:   HANDLE_UNARY_CONST_FUNC(log); break;
            case cAcosh: HANDLE_UNARY_CONST_FUNC(fp_acosh); break;
            case cAsinh: HANDLE_UNARY_CONST_FUNC(fp_asinh); break;
            case cAtanh: HANDLE_UNARY_CONST_FUNC(fp_atanh); break;
            case cAcos: HANDLE_UNARY_CONST_FUNC(acos); break;
            case cAsin: HANDLE_UNARY_CONST_FUNC(asin); break;
            case cAtan: HANDLE_UNARY_CONST_FUNC(atan); break;
            case cCosh: HANDLE_UNARY_CONST_FUNC(cosh); break;
            case cSinh: HANDLE_UNARY_CONST_FUNC(sinh); break;
            case cTanh: HANDLE_UNARY_CONST_FUNC(tanh); break;
            case cSin: HANDLE_UNARY_CONST_FUNC(sin); break;
            case cCos: HANDLE_UNARY_CONST_FUNC(cos); break;
            case cTan: HANDLE_UNARY_CONST_FUNC(tan); break;
            case cCeil: HANDLE_UNARY_CONST_FUNC(ceil); break;
            case cFloor: HANDLE_UNARY_CONST_FUNC(floor); break;
            case cInt:
                if(Params[0].param->IsImmed())
                    { const_value = floor(Params[0].param->GetImmed() + 0.5);
                      goto ReplaceTreeWithConstValue; }
                break;
            case cLog2:
                if(Params[0].param->IsImmed())
                    { const_value = log(Params[0].param->GetImmed()) * CONSTANT_L2I;
                      goto ReplaceTreeWithConstValue; }
                break;
            case cLog10:
                if(Params[0].param->IsImmed())
                    { const_value = log(Params[0].param->GetImmed()) * CONSTANT_L10I;
                      goto ReplaceTreeWithConstValue; }
                break;

            case cAtan2:
            {
                /* Range based optimizations for (y,x):
                 * If y is +0 and x <= -0, +pi is returned
                 * If y is -0 and x <= -0, -pi is returned (assumed never happening)
                 * If y is +0 and x >= +0, +0 is returned
                 * If y is -0 and x >= +0, -0 is returned  (assumed never happening)
                 * If x is +-0 and y < 0, -pi/2 is returned
                 * If x is +-0 and y > 0, +pi/2 is returned
                 * Otherwise, perform constant folding when available
                 * If we know x <> 0, convert into atan(y / x)
                 *   TODO: Figure out whether the above step is wise
                 *         It allows e.g. atan2(6*x, 3*y) -> atan(2*x/y)
                 *         when we know y != 0
                 */
                MinMaxTree p0 = Params[0].param->CalculateResultBoundaries();
                MinMaxTree p1 = Params[1].param->CalculateResultBoundaries();
                if(p0.has_min && p0.has_max && p0.min == 0.0)
                {
                    if(p1.has_max && p1.max < 0)
                        { const_value = CONSTANT_PI; goto ReplaceTreeWithConstValue; }
                    if(p1.has_max && p1.max >= 0.0)
                        { const_value = p0.min; goto ReplaceTreeWithConstValue; }
                }
                if(p1.has_min && p1.has_max && p1.min == 0.0)
                {
                    if(p0.has_max && p0.max < 0)
                        { const_value = -CONSTANT_PIHALF; goto ReplaceTreeWithConstValue; }
                    if(p0.has_min && p0.min > 0)
                        { const_value =  CONSTANT_PIHALF; goto ReplaceTreeWithConstValue; }
                }
                if(Params[0].param->IsImmed()
                && Params[1].param->IsImmed())
                    { const_value = atan2(Params[0].param->GetImmed(),
                                          Params[1].param->GetImmed());
                      goto ReplaceTreeWithConstValue; }
              #if 0
                if((p1.has_min && p1.min > 0.0)
                || (p1.has_max && p1.max < NEGATIVE_MAXIMUM))
                {
                    // Convert into a division
                    CodeTreeP subtree = new CodeTree;
                    Params[1].sign = true; /* FIXME: Not appropriate anymore */
                    for(size_t a=0; a<Params.size(); ++a)
                        Params[a].param->Parent = &*subtree;
                    subtree->Opcode = cMul;
                    subtree->Params.swap(Params); // subtree = y/x
                    subtree->ConstantFolding();
                    subtree->Sort();
                    subtree->Rehash(false);
                    Opcode = cAtan;
                    AddParam(Param(subtree, false)); // we = atan(y/x)
                }
              #endif
                break;
            }

            case cPow:
            {
                if(Params[0].param->IsImmed()
                && Params[1].param->IsImmed())
                    { const_value = pow(Params[0].param->GetImmed(),
                                        Params[1].param->GetImmed());
                      goto ReplaceTreeWithConstValue; }
                if(Params[1].param->IsImmed()
                && Params[1].param->GetImmed() == 1.0)
                {
                    // x^1 = x
                    goto ReplaceTreeWithParam0;
                }
                if(Params[0].param->IsImmed()
                && Params[0].param->GetImmed() == 1.0)
                {
                    // 1^x = 1
                    goto ReplaceTreeWithOne;
                }
                break;
            }

            case cMod:
            {
                /* Can more be done than this? */
                if(Params[0].param->IsImmed()
                && Params[1].param->IsImmed())
                    { const_value = fmod(Params[0].param->GetImmed(),
                                         Params[1].param->GetImmed());
                      goto ReplaceTreeWithConstValue; }
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.cc and thus
             * they will never occur in the calling context:
             */
            case cNeg: // converted into cAdd ~x
            case cInv: // converted into cMul ~x
            case cDiv: // converted into cMul ~x
            case cRDiv: // similar to above
            case cSub: // converted into cAdd ~x
            case cRSub: // similar to above
            case cRad: // converted into cMul x CONSTANT_RD
            case cDeg: // converted into cMul x CONSTANT_DR
            case cSqr: // converted into cMul x x
            case cExp: // converted into cPow CONSTANT_E x
            case cExp2: // converted into cPow 2.0 x
            case cSqrt: // converted into cPow x 0.5
            case cRSqrt: // converted into cPow x -0.5
            case cCot: // converted into cMul ~(cTan x)
            case cSec: // converted into cMul ~(cCos x)
            case cCsc: // converted into cMul ~(cSin x)
                break; /* Should never occur */

            /* Opcodes that do not occur in the tree for other reasons */
            case cDup:
            case cFetch:
            case cPopNMov:
            case cNop:
            case cJump:
            case VarBegin:
                break; /* Should never occur */
            /* Opcodes that we can't do anything about */
            case cPCall:
            case cFCall:
#         ifndef FP_DISABLE_EVAL
            case cEval:
#endif
                break;
        }
        /*
        if(Parent)
            Parent->ConstantFolding();

        */
    }
}

namespace FPoptimizer_Grammar
{
    static double GetPackConst(size_t index)
    {
        double res = pack.clist[index];
    #if 0
        if(res == FPOPT_NAN_CONST)
        {
        #ifdef NAN
            return NAN;
        #else
            return 0.0; // Should be 0.0/0.0, but some compilers don't like that
        #endif
        }
    #endif
        return res;
    }

    /* A helper for std::equal_range */
    struct OpcodeRuleCompare
    {
        bool operator() (const FPoptimizer_CodeTree::CodeTree& tree, const Rule& rule) const
        {
            /* If this function returns true, len=half.
             */

            if(tree.Opcode != rule.func.opcode)
                return tree.Opcode < rule.func.opcode;

            if(tree.Params.size() < rule.n_minimum_params)
            {
                // Tree has fewer params than required?
                return true; // Failure
            }
            return false;
        }
        bool operator() (const Rule& rule, const FPoptimizer_CodeTree::CodeTree& tree) const
        {
            /* If this function returns true, rule will be excluded from the equal_range
             */

            if(rule.func.opcode != tree.Opcode)
                return rule.func.opcode < tree.Opcode;

            if(rule.n_minimum_params < tree.Params.size())
            {
                // Tree has more params than the pattern has?
                switch(pack.mlist[rule.func.index].type)
                {
                    case PositionalParams:
                    case SelectedParams:
                        return true; // Failure
                    case AnyParams:
                        return false; // Not a failure
                }
            }
            return false;
        }
    };

#ifdef DEBUG_SUBSTITUTIONS
    static const char ImmedHolderNames[4][2]  = {"%","&"};
    static const char NamedHolderNames[10][2] = {"x","y","z","a","b","c","d","e","f","g"};
#endif

    /* Apply the grammar to a given CodeTree */
    bool Grammar::ApplyTo(
        FPoptimizer_CodeTree::CodeTree& tree,
        bool recursion) const
    {
        bool changed = false;

        //        recursion=recursion;

        if(tree.OptimizedUsing != this)
        {
            /* First optimize all children */
            tree.ConstantFolding();

            for(size_t a=0; a<tree.Params.size(); ++a)
            {
                if( ApplyTo( *tree.Params[a].param, true ) )
                {
                    changed = true;
                }
            }

            if(changed)
            {
                // Give the parent node a rerun at optimization
                return true;
            }

            /* Figure out which rules _may_ match this tree */
            typedef const Rule* ruleit;

            std::pair<ruleit, ruleit> range
                = MyEqualRange(pack.rlist + this->index,
                               pack.rlist + this->index + this->count,
                               tree,
                               OpcodeRuleCompare());

#ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Input (Grammar #"
                      << (this - pack.glist)
                      << ", " << FP_GetOpcodeName(tree.Opcode)
                      << "[" << tree.Params.size()
                      << "], rules "
                      << (range.first - pack.rlist)
                      << ".."
                      << (range.second - pack.rlist)
                      << ": ";
            DumpTree(tree);
            std::cout << "\n" << std::flush;
#endif

            while(range.first != range.second)
            {
                /* Check if this rule matches */
                if(range.first->ApplyTo(tree))
                {
                    changed = true;
                    break;
                }
                ++range.first;
            }

#ifdef DEBUG_SUBSTITUTIONS
            std::cout << (changed ? "Changed." : "No changes.");
            std::cout << "\n" << std::flush;
#endif

            if(!changed)
            {
                tree.OptimizedUsing = this;
            }
        }
        else
        {
#ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Already optimized:  ";
            DumpTree(tree);
            std::cout << "\n" << std::flush;
#endif
        }

#ifdef DEBUG_SUBSTITUTIONS
        if(!recursion)
        {
            std::cout << "Output: ";
            DumpTree(tree);
            std::cout << "\n" << std::flush;
        }
#endif
        return changed;
    }

    /* Store information about a potential match,
     * in order to iterate through candidates
     */
    struct MatchedParams::CodeTreeMatch
    {
        // Which parameters were matched -- these will be replaced if AnyParams are used
        std::vector<size_t> param_numbers;

        // Which values were saved for ImmedHolders?
        std::map<unsigned, double> ImmedMap;
        // Which codetrees were saved for each NameHolder? And how many?
            struct NamedItem
            {
                fphash_t hash;
                size_t   howmany;
                size_t   n_synthesized;

                NamedItem() { howmany = 0; n_synthesized = 0; }
                NamedItem(fphash_t h,size_t m) { hash = h; howmany = m; n_synthesized = 0; }
            };
        std::map<unsigned, NamedItem> NamedMap;
        // Which codetrees were saved for each RestHolder?
        std::map<unsigned,
          std::vector<fphash_t> > RestMap;

        // Examples of each codetree
        std::map<fphash_t, FPoptimizer_CodeTree::CodeTreeP> trees;

        CodeTreeMatch() : param_numbers(), ImmedMap(), NamedMap(), RestMap() { }
    };

#ifdef DEBUG_SUBSTITUTIONS
    void DumpMatch(const Function& input,
                   const FPoptimizer_CodeTree::CodeTree& tree,
                   const MatchedParams& replacement,
                   const MatchedParams::CodeTreeMatch& matchrec,
                   bool DidMatch=true);
    void DumpFunction(const Function& input);
    void DumpParam(const ParamSpec& p);
    void DumpParams(const MatchedParams& mitem);
#endif

    /* Apply the rule to a given CodeTree */
    bool Rule::ApplyTo(
        FPoptimizer_CodeTree::CodeTree& tree) const
    {
        const Function&      input  = func;
        const MatchedParams& repl   = pack.mlist[repl_index];

        if(input.opcode == tree.Opcode)
        {
            for(unsigned long match_index=0; ; ++match_index)
            {
                MatchedParams::CodeTreeMatch matchrec;
                MatchResultType mr =
                    pack.mlist[input.index].Match(tree, matchrec,match_index, false);
                if(!mr.found && mr.has_more) continue;
                if(!mr.found) break;

    #ifdef DEBUG_SUBSTITUTIONS
                DumpMatch(input, tree, repl, matchrec);
    #endif

                const MatchedParams& params = pack.mlist[input.index];
                switch(type)
                {
                    case ReplaceParams:
                        repl.ReplaceParams(tree, params, matchrec);
    #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "  ParmReplace: ";
                        DumpTree(tree);
                        std::cout << "\n" << std::flush;
                        DumpHashes(tree);
    #endif
                        return true;
                    case ProduceNewTree:
                        repl.ReplaceTree(tree,   params, matchrec);
    #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "  TreeReplace: ";
                        DumpTree(tree);
                        std::cout << "\n" << std::flush;
                        DumpHashes(tree);
    #endif
                        return true;
                }
                break; // should be unreachable
            }
        }
        #ifdef DEBUG_SUBSTITUTIONS
        // Report mismatch
        MatchedParams::CodeTreeMatch matchrec;
        DumpMatch(input, tree, repl, matchrec, false);
        #endif
        return false;
    }


    /* Match the given function to the given CodeTree.
     */
    MatchResultType Function::Match(
        FPoptimizer_CodeTree::CodeTree& tree,
        MatchedParams::CodeTreeMatch& match,
        unsigned long match_index) const
    {
        if(opcode != tree.Opcode) return NoMatch;
        return pack.mlist[index].Match(tree, match, match_index, true);
    }


    /* This struct is used by MatchedParams::Match() for backtracking. */
    struct ParamMatchSnapshot
    {
        MatchedParams::CodeTreeMatch snapshot;
                                    // Snapshot of the state so far
        size_t            parampos; // Which position was last chosen?
        std::vector<bool> used;     // Which params were allocated?

        size_t            matchpos;
    };

    /* Match the given list of ParamSpecs using the given ParamMatchingType
     * to the given CodeTree.
     * The CodeTree is already assumed to be a function type
     * -- i.e. it is assumed that the caller has tested the Opcode of the tree.
     */
    MatchResultType MatchedParams::Match(
        FPoptimizer_CodeTree::CodeTree& tree,
        MatchedParams::CodeTreeMatch& match,
        unsigned long match_index,
        bool recursion) const
    {
        /*        match_index is a feature for backtracking.
         *
         *        For example,
         *          cMul (cAdd x) (cAdd x)
         *        Applied to:
         *          (a+b)*(c+b)
         *
         *        Match (cAdd x) to (a+b) may first capture "a" into "x",
         *        and then Match(cAdd x) for (c+b) will fail,
         *        because there's no "a" there.
         *
         *        However, match_index can be used to indicate that the
         *        _second_ matching will be used, so that "b" will be
         *        captured into "x".
         */


        /* First, check if the tree has any chances of matching... */
        /* Figure out what we need. */
        struct Needs
        {
            struct Needs_Pol
            {
                int SubTrees; // This many subtrees
                int Others;   // This many others (namedholder)
                unsigned SubTreesDetail[VarBegin]; // This many subtrees of each opcode type

                Needs_Pol(): SubTrees(0), Others(0), SubTreesDetail()
                {
                }
            } polarity[2]; // 0=positive, 1=negative

            int Immeds;      // This many immeds

            Needs(): polarity(), Immeds() { }
        } NeedList;

        // Figure out what we need
        size_t minimum_need = 0;
        for(unsigned a=0; a<count; ++a)
        {
            const ParamSpec& param = pack.plist[index+a];
            Needs::Needs_Pol& needs = NeedList.polarity[param.sign];
            switch(param.opcode)
            {
                case SubFunction:
                    needs.SubTrees += 1;
                    assert( pack.flist[param.index].opcode < VarBegin );
                    needs.SubTreesDetail[ pack.flist[param.index].opcode ] += 1;
                    ++minimum_need;
                    break;
                case NumConstant:
                case ImmedHolder:
                default: // GroupFunction:
                    NeedList.Immeds += 1;
                    ++minimum_need;
                    break;
                case NamedHolder:
                    needs.Others += param.minrepeat;
                    ++minimum_need;
                    break;
                case RestHolder:
                    break;
            }
        }
        if(tree.Params.size() < minimum_need)
        {
            // Impossible to satisfy
            return NoMatch;
        }

        // Figure out what we have (note: we already assume that the opcode of the tree matches!)
        for(size_t a=0; a<tree.Params.size(); ++a)
        {
            Needs::Needs_Pol& needs = NeedList.polarity[tree.Params[a].sign];
            unsigned opcode = tree.Params[a].param->Opcode;
            switch(opcode)
            {
                case cImmed:
                    if(NeedList.Immeds > 0) NeedList.Immeds -= 1;
                    else needs.Others -= 1;
                    break;
                case cVar:
                case cFCall:
                case cPCall:
                    needs.Others -= 1;
                    break;
                default:
                    assert( opcode < VarBegin );
                    if(needs.SubTrees > 0
                    && needs.SubTreesDetail[opcode] > 0)
                    {
                        needs.SubTrees -= 1;
                        needs.SubTreesDetail[opcode] -= 1;
                    }
                    else needs.Others -= 1;
            }
        }

        // Check whether all needs were satisfied
        if(NeedList.Immeds > 0
        || NeedList.polarity[0].SubTrees > 0
        || NeedList.polarity[0].Others > 0
        || NeedList.polarity[1].SubTrees > 0
        || NeedList.polarity[1].Others > 0)
        {
            // Something came short, impossible to satisfy.
            return NoMatch;
        }

        if(type != AnyParams)
        {
            if(NeedList.Immeds < 0
            || NeedList.polarity[0].SubTrees < 0
            || NeedList.polarity[0].Others < 0
            || NeedList.polarity[1].SubTrees < 0
            || NeedList.polarity[1].Others < 0
            || count != tree.Params.size())
            {
                // Something was too much.
                return NoMatch;
            }
        }

        TransformationType transf = None;
        switch(tree.Opcode)
        {
            case cAdd: transf = Negate; break;
            case cMul: transf = Invert; break;
            case cAnd:
            case cOr:  transf = NotThe; break;
        }

        switch(type)
        {
            case PositionalParams:
            {
                /*DumpTree(tree);
                std::cout << "<->";
                DumpParams(*this);
                std::cout << " -- ";*/

                std::vector<MatchPositionSpec<CodeTreeMatch> > specs;
                specs.reserve(count);
                //fprintf(stderr, "Enter loop %lu\n", match_index);
                for(unsigned a=0; a<count; ++a)
                {
                    specs.resize(a+1);

                PositionalParamsMatchingLoop:;
                    // Match this parameter.
                    MatchResultType mr = pack.plist[index+a].Match(
                        *tree.Params[a].param, match,
                        tree.Params[a].sign ? transf : None,
                        specs[a].roundno);

                    specs[a].done = !mr.has_more;

                    // If it was not found, backtrack...
                    if(!mr.found)
                    {
                    LoopThisRound:
                        while(specs[a].done)
                        {
                            // Backtrack
                            if(a <= 0) return NoMatch; //
                            specs.resize(a);
                            --a;
                            match = specs[a].data;
                        }
                        ++specs[a].roundno;
                        goto PositionalParamsMatchingLoop;
                    }
                    // If found...
                    if(!recursion)
                        match.param_numbers.push_back(a);
                    specs[a].data = match;

                    if(a == count-1U && match_index > 0)
                    {
                        // Skip this match
                        --match_index;
                        goto LoopThisRound;
                    }
                }
                /*std::cout << " yay?\n";*/
                // Match = no mismatch.
                bool final_try = true;
                for(unsigned a=0; a<count; ++a)
                    if(!specs[a].done) { final_try = false; break; }
                //fprintf(stderr, "Exit  loop %lu\n", match_index);
                return MatchResultType(true, !final_try);
            }
            case AnyParams:
            case SelectedParams:
            {
                const size_t n_tree_params = tree.Params.size();

                unsigned N_PositiveRestHolders = 0;
                unsigned N_NegativeRestHolders = 0;
                for(unsigned a=0; a<count; ++a)
                {
                    const ParamSpec& param = pack.plist[index+a];
                    if(param.opcode == RestHolder)
                    {
                        if(param.sign)
                            ++N_NegativeRestHolders;
                        else
                            ++N_PositiveRestHolders;
                    }
                }

                bool HasRestHolders = N_PositiveRestHolders || N_NegativeRestHolders;

                #ifdef DEBUG_SUBSTITUTIONS
                if((type == AnyParams) && recursion && !HasRestHolders)
                {
                    std::cout << "Recursed AnyParams with no RestHolders?\n";
                    DumpParams(*this);
                }
                #endif

                if(!HasRestHolders && recursion && count != n_tree_params)
                {
                    /*DumpTree(tree);
                    std::cout << "<->";
                    DumpParams(*this);
                    std::cout << " -- fail due to recursion&&count!=n_tree_params";*/
                    return NoMatch; // Impossible match.
                }

                /*std::cout << "Matching ";
                DumpTree(tree); std::cout << " with ";
                DumpParams(*this);
                std::cout << " , match_index=" << match_index << "\n" << std::flush;*/

                std::vector<ParamMatchSnapshot> position(count);
                std::vector<bool>               used(n_tree_params);

                unsigned p=0;

                for(; p<count; ++p)
                {
                    position[p].snapshot  = match;
                    position[p].parampos  = 0;
                    position[p].matchpos  = 0;
                    position[p].used      = used;

                    //fprintf(stderr, "posA: p=%u count=%u\n", p, count);

                backtrack:
                  {
                    if(pack.plist[index+p].opcode == RestHolder)
                    {
                        // RestHolders always match. They're filled afterwards.
                        position[p].parampos = n_tree_params;
                        position[p].matchpos = 0;
                        continue;
                    }

                    size_t whichparam = position[p].parampos;
                    size_t whichmatch = position[p].matchpos;

                    /* a          = param index in the syntax specification
                     * whichparam = param index in the tree received from parser
                     */

                    /*fprintf(stderr, "posB: p=%u, whichparam=%lu, whichmatch=%lu\n",
                        p,whichparam,whichmatch);*/
                    while(whichparam < n_tree_params)
                    {
                        if(used[whichparam])
                        {
                        NextParamNumber:
                            ++whichparam;
                            whichmatch = 0;
                            continue;
                        NextMatchNumber:
                            ++whichmatch;
                        }

                        /*std::cout << "Maybe [" << p << "]:";
                        DumpParam(pack.plist[index+p]);
                        std::cout << " <-> ";
                        if(tree.Params[whichparam].sign) std::cout << '~';
                        DumpTree(*tree.Params[whichparam].param);
                        std::cout << "...?\n" << std::flush;*/

                        MatchResultType mr = pack.plist[index+p].Match(
                            *tree.Params[whichparam].param, match,
                            tree.Params[whichparam].sign ? transf : None,
                            whichmatch);

                        /*std::cout << "In ";
                        DumpTree(tree); std::cout << std::flush;
                        fprintf(stderr, ", trying param %lu, match %lu (matchindex %lu); got %s,%s: ",
                            whichparam,whichmatch, match_index,
                            mr.found?"found":"not found",
                            mr.has_more?"more":"no more"); fflush(stderr);
                        DumpParam(pack.plist[index+p]); std::cout << "\n" << std::flush;*/

                        if(!mr.found)
                        {
                        NextParamTest:
                            if(!mr.has_more) goto NextParamNumber;
                            goto NextMatchNumber;
                        }

                        /*std::cout << "woo... " << a << ", " << b << "\n";*/
                        /* NamedHolders require a special treatment,
                         * because a repetition count may be issued
                         * for them.
                         */
                        if(pack.plist[index+p].opcode == NamedHolder)
                        {
                            // Verify the MinRepeat & AnyRepeat case
                            unsigned MinRepeat = pack.plist[index+p].minrepeat;
                            bool AnyRepeat     = pack.plist[index+p].anyrepeat;
                            unsigned HadRepeat = 1;

                            for(size_t repeat_pos = whichparam+1;
                                repeat_pos < n_tree_params && (HadRepeat < MinRepeat || AnyRepeat);
                                ++repeat_pos)
                            {
                                /*fprintf(stderr, "Req @ %lu = %d:%16lX, got @ %lu = %d:%16lX\n",
                                    whichparam, tree.Params[whichparam].sign,
                                                tree.Params[whichparam].param->Hash,
                                    repeat_pos, tree.Params[repeat_pos].sign,
                                                tree.Params[repeat_pos].param->Hash);*/

                                if(tree.Params[repeat_pos].IsIdenticalTo(
                                   tree.Params[whichparam])
                                && !used[repeat_pos])
                                {
                                    ++HadRepeat;
                                }
                            }
                            /*fprintf(stderr, "Got repeat %u, needs %u\n", HadRepeat,MinRepeat);*/
                            if(HadRepeat < MinRepeat)
                            {
                                match = position[p].snapshot;
                                used  = position[p].used;
                                goto NextParamTest; // No sufficient repeat count here
                            }

                            used[whichparam] = true;
                            if(!recursion) match.param_numbers.push_back(whichparam);

                            HadRepeat = 1;
                            for(size_t repeat_pos = whichparam+1;
                                repeat_pos < n_tree_params && (HadRepeat < MinRepeat || AnyRepeat);
                                ++repeat_pos)
                            {
                                if(tree.Params[repeat_pos].IsIdenticalTo(
                                   tree.Params[whichparam])
                                && !used[repeat_pos])
                                {
                                    ++HadRepeat;
                                    used[repeat_pos] = true;
                                    if(!recursion) match.param_numbers.push_back(repeat_pos);
                                }
                            }
                            if(AnyRepeat)
                                match.NamedMap[pack.plist[index+p].index].howmany = HadRepeat;
                        }
                        else
                        {
                            used[whichparam] = true;
                            if(!recursion) match.param_numbers.push_back(whichparam);
                        }
                        position[p].parampos = mr.has_more ? whichparam : (whichparam+1);
                        position[p].matchpos = mr.has_more ? (whichmatch+1) : 0;
                        goto ok;
                    }

                    /*DumpParam(param);
                    std::cout << " didn't match anything in ";
                    DumpTree(tree);
                    std::cout << "\n";*/
                  }

                    // No match for this param, try backtracking.
                DiscardedThisAttempt:
                    while(p > 0)
                    {
                        --p;
                        ParamMatchSnapshot& prevpos = position[p];
                        if(prevpos.parampos < n_tree_params)
                        {
                            // Try another combination.
                            match = prevpos.snapshot;
                            used  = prevpos.used;
                            goto backtrack;
                        }
                    }
                    // If we cannot backtrack, break. No possible match.
                    /*if(!recursion)
                        std::cout << "Drats!\n";*/
                    if(match_index == 0)
                        return NoMatch;
                    break;
                ok:;
                    /*if(!recursion)
                        std::cout << "Match for param " << a << " at " << b << std::endl;*/

                    if(p == count-1U && match_index > 0)
                    {
                        // Skip this match
                        --match_index;
                        goto DiscardedThisAttempt;
                    }
                }
                /*fprintf(stderr, "End loop, match_index=%lu\n", match_index); fflush(stderr);*/

                /* We got a match. */

                // If the rule cares about the balance of
                // negative restholdings versus positive restholdings,
                // verify them.
                if(balance != BalanceDontCare)
                {
                    unsigned n_pos_restholdings = 0;
                    unsigned n_neg_restholdings = 0;

                    for(unsigned a=0; a<count; ++a)
                    {
                        const ParamSpec& param = pack.plist[index+a];
                        if(param.opcode == RestHolder)
                        {
                            for(size_t b=0; b<n_tree_params; ++b)
                                if(tree.Params[b].sign == param.sign && !used[b])
                                {
                                    if(param.sign)
                                        n_neg_restholdings += 1;
                                    else
                                        n_pos_restholdings += 1;
                                }
                        }
                    }
                    switch(balance)
                    {
                        case BalanceMoreNeg:
                            if(n_neg_restholdings <= n_pos_restholdings) return NoMatch;
                            break;
                        case BalanceMorePos:
                            if(n_pos_restholdings <= n_neg_restholdings) return NoMatch;
                            break;
                        case BalanceEqual:
                            if(n_pos_restholdings != n_neg_restholdings) return NoMatch;
                            break;
                        case BalanceDontCare: ;
                    }
                }

                unsigned pos_rest_remain = N_PositiveRestHolders;
                unsigned neg_rest_remain = N_NegativeRestHolders;

                // Verify if we have RestHolder constraints.
                for(unsigned a=0; a<count; ++a)
                {
                    const ParamSpec& param = pack.plist[index+a];
                    if(param.opcode == RestHolder)
                    {
                        std::map<unsigned, std::vector<fphash_t> >::iterator
                            i = match.RestMap.lower_bound(param.index);

                        if(i != match.RestMap.end() && i->first == param.index)
                        {
                            unsigned& n_remaining_restholders_of_this_kind =
                                param.sign ? neg_rest_remain : pos_rest_remain;
                            /*fprintf(stderr, "Does restholder %u match in", param.index);
                            fflush(stderr); DumpTree(tree); std::cout << "? " << std::flush;*/

                            const std::vector<fphash_t>& RefRestList = i->second;
                            for(size_t r=0; r<RefRestList.size(); ++r)
                            {
                                for(size_t b=0; b<n_tree_params; ++b)
                                    if(tree.Params[b].sign == param.sign
                                    && !used[b]
                                    && tree.Params[b].param->Hash == RefRestList[r]
                                    && tree.Params[b].param->IsIdenticalTo(
                                        * match.trees.find(RefRestList[r])->second )
                                      )
                                    {
                                        used[b] = true;
                                        goto SatisfiedRestHolder;
                                    }
                                // Unsatisfied RestHolder constraint
                                /*fprintf(stderr, "- no\n");*/
                                p=count-1;
                                goto DiscardedThisAttempt;
                            SatisfiedRestHolder:;
                            }
                            --n_remaining_restholders_of_this_kind;
                            /*fprintf(stderr, "- yes\n");*/
                        }
                    }
                }

                // Now feed any possible RestHolders the remaining parameters.
                bool more_restholder_options = false;
                for(unsigned a=0; a<count; ++a)
                {
                    const ParamSpec& param = pack.plist[index+a];
                    if(param.opcode == RestHolder)
                    {
                        std::map<unsigned, std::vector<fphash_t> >::iterator
                            i = match.RestMap.lower_bound(param.index);
                        if(i != match.RestMap.end() && i->first == param.index) continue;

                        std::vector<fphash_t>& RestList = match.RestMap[param.index]; // mark it up

                        unsigned& n_remaining_restholders_of_this_kind =
                            param.sign ? neg_rest_remain : pos_rest_remain;

                        unsigned n_remaining_params = 0;
                        for(size_t b=0; b<n_tree_params; ++b)
                            if(tree.Params[b].sign == param.sign && !used[b])
                                ++n_remaining_params;

                        /*fprintf(stderr, "[index %lu] For restholder %u, %u remains, %u remaining of kind\n",
                            match_index,
                            (unsigned)param.index, (unsigned)n_remaining_params,
                            (unsigned)n_remaining_restholders_of_this_kind);
                            fflush(stderr);*/

                        if(n_remaining_params > 0)
                        {
                            if(n_remaining_params > 8) n_remaining_params = 8;
                            unsigned n_remaining_combinations = 1 << n_remaining_params;

                            unsigned n_options = n_remaining_restholders_of_this_kind > 1
                                ? n_remaining_combinations
                                : 1;
                            size_t selection = n_remaining_combinations - 1;
                            if(n_options > 1)
                            {
                                --n_options;
                                selection = match_index % (n_options); ++selection;
                                match_index /= n_options;
                            }
                            if(selection+1 < n_options) more_restholder_options = true;

                            /*fprintf(stderr, "- selected %u/%u\n", selection, n_options); fflush(stderr);*/

                            unsigned matchbit = 1;
                            for(size_t b=0; b<n_tree_params; ++b)
                                if(tree.Params[b].sign == param.sign && !used[b])
                                {
                                    if(selection & matchbit)
                                    {
                                        /*fprintf(stderr, "- uses param %lu\n", b);*/
                                        if(!recursion)
                                            match.param_numbers.push_back(b);
                                        fphash_t hash = tree.Params[b].param->Hash;
                                        RestList.push_back(hash);
                                        match.trees.insert(
                                            std::make_pair(hash, tree.Params[b].param) );

                                        used[b] = true;
                                    }
                                    if(matchbit < 0x80U) matchbit <<= 1;
                                }
                        }
                        --n_remaining_restholders_of_this_kind;
                    }
                }
                /*std::cout << "Returning match for ";
                DumpTree(tree);
                std::cout << "\n               with ";
                DumpParams(*this); std::cout << std::flush;
                fprintf(stderr, ", %s hope for more (now %lu)\n",
                    more_restholder_options ? "with" : "without", match_index); fflush(stderr);*/
                return more_restholder_options ? FoundSomeMatch : FoundLastMatch;
            }
        }
        return NoMatch;
    }

    MatchResultType ParamSpec::Match(
        FPoptimizer_CodeTree::CodeTree& tree,
        MatchedParams::CodeTreeMatch& match,
        TransformationType transf,
        unsigned long match_index) const
    {
        assert(opcode != RestHolder); // RestHolders are supposed to be handled by the caller

        switch(OpcodeType(opcode))
        {
            case NumConstant:
            {
                if(!tree.IsImmed()) return NoMatch;
                double res = tree.GetImmed();
                if(transformation == Negate) res = -res;
                if(transformation == Invert) res = 1/res;
                double res2 = GetPackConst(index);
                if(transf == Negate) res2 = -res2;
                if(transf == Invert) res2 = 1/res2;
                if(transf == NotThe) res2 = res2 != 0;
                /*std::cout << std::flush;
                fprintf(stderr, "Comparing %.20f and %.20f\n", res, res2);
                fflush(stderr);*/
                if(!FloatEqual(res, res2)) return NoMatch;
                return FoundLastMatch; // Previously unknown NumConstant, good
            }
            case ImmedHolder:
            {
                if(!tree.IsImmed()) return NoMatch;
                double res = tree.GetImmed();

                switch( ImmedConstraint_Value(count & ValueMask) )
                {
                    case ValueMask: break;
                    case Value_AnyNum: break;
                    case Value_EvenInt:
                        if(!FloatEqual(res, (double)(long)(res))) return NoMatch;
                        if( (long)(res) % 2 != 0) return NoMatch;
                        break;
                    case Value_OddInt:
                        if(!FloatEqual(res, (double)(long)(res))) return NoMatch;
                        if( (long)(res) % 2 == 0) return NoMatch;
                        break;
                    case Value_IsInteger:
                        if(!FloatEqual(res, (double)(long)(res))) return NoMatch;
                        break;
                    case Value_NonInteger:
                        if(FloatEqual(res, (double)(long)(res))) return NoMatch;
                        break;
                }
                switch( ImmedConstraint_Sign(count & SignMask) )
                {
                    /*case SignMask: break;*/
                    case Sign_AnySign: break;
                    case Sign_Positive:
                        if(res < 0.0)  return NoMatch;
                        break;
                    case Sign_Negative:
                        if(res >= 0.0) return NoMatch;
                        break;
                    case Sign_NoIdea:
                        return NoMatch;
                }
                switch( ImmedConstraint_Oneness(count & OnenessMask) )
                {
                    case OnenessMask: break;
                    case Oneness_Any: break;
                    case Oneness_One:
                        if(!FloatEqual(fabs(res), 1.0)) return NoMatch;
                        break;
                    case Oneness_NotOne:
                        if(FloatEqual(fabs(res), 1.0)) return NoMatch;
                        break;
                }

                if(transformation == Negate) res = -res;
                if(transformation == Invert) res = 1/res;
                std::map<unsigned, double>::iterator
                    i = match.ImmedMap.lower_bound(index);
                if(i != match.ImmedMap.end() && i->first == index)
                {
                    double res2 = i->second;
                    if(transf == Negate) res2 = -res2;
                    if(transf == Invert) res2 = 1/res2;
                    if(transf == NotThe) res2 = res2 != 0;
                    /*std::cout << std::flush;
                    fprintf(stderr, "Comparing %.20f and %.20f\n", res, res2);
                    fflush(stderr);*/
                    return FloatEqual(res, res2) ? FoundLastMatch : NoMatch;
                }
                if(sign != (transf != None)) return NoMatch;

                unsigned ind = (unsigned)index;
                match.ImmedMap.insert(i, std::make_pair(ind, res));
                return FoundLastMatch; // Previously unknown ImmedHolder, good
            }
            case NamedHolder:
            {
                if(sign != (transf != None)) return NoMatch;
                std::map<unsigned, MatchedParams::CodeTreeMatch::NamedItem>::iterator
                    i = match.NamedMap.lower_bound(index);
                if(i != match.NamedMap.end() && i->first == index)
                {
                    /*fprintf(stderr, "NamedHolder found: %16lX -- tested against %16lX\n", i->second.first, tree.Hash);*/
                    if(tree.Hash == i->second.hash
                    && tree.IsIdenticalTo(* match.trees.find(i->second.hash)->second)
                      )
                        return FoundLastMatch;
                    else
                        return NoMatch;
                }

                switch( ImmedConstraint_Value(count & ValueMask) )
                {
                    case ValueMask: break;
                    case Value_AnyNum: break;
                    case Value_EvenInt:
                        if(!tree.IsAlwaysParity(false)) return NoMatch;
                        break;
                    case Value_OddInt:
                        if(!tree.IsAlwaysParity(true)) return NoMatch;
                        break;
                    case Value_IsInteger:
                        if(!tree.IsAlwaysInteger()) return NoMatch;
                        break;
                    case Value_NonInteger:
                        if(tree.IsAlwaysInteger()) return NoMatch;
                        break;
                }
                switch( ImmedConstraint_Sign(count & SignMask) )
                {
                    /*case SignMask: break;*/
                    case Sign_AnySign: break;
                    case Sign_Positive:
                        if(!tree.IsAlwaysSigned(true)) return NoMatch;
                        break;
                    case Sign_Negative:
                        if(!tree.IsAlwaysSigned(false)) return NoMatch;
                        break;
                    case Sign_NoIdea:
                        if(tree.IsAlwaysSigned(false)) return NoMatch;
                        if(tree.IsAlwaysSigned(true)) return NoMatch;
                        break;
                }
                switch( ImmedConstraint_Oneness(count & OnenessMask) )
                {
                    case OnenessMask: break;
                    case Oneness_Any: break;
                    case Oneness_One:    return NoMatch;
                    case Oneness_NotOne: return NoMatch;
                }

                match.NamedMap.insert(i,
                    std::make_pair(index,
                        MatchedParams::CodeTreeMatch::NamedItem(tree.Hash,1) ));
                match.trees.insert(std::make_pair(tree.Hash, &tree));
                return FoundLastMatch; // Previously unknown NamedHolder, good
            }
            case RestHolder:
            {
                break;
            }
            case SubFunction:
            {
                if(sign != (transf != None)) return NoMatch;

                switch( ImmedConstraint_Value(count & ValueMask) )
                {
                    case ValueMask: break;
                    case Value_AnyNum: break;
                    case Value_EvenInt:
                        if(!tree.IsAlwaysParity(false)) return NoMatch;
                        break;
                    case Value_OddInt:
                        if(!tree.IsAlwaysParity(true)) return NoMatch;
                        break;
                    case Value_IsInteger:
                        if(!tree.IsAlwaysInteger()) return NoMatch;
                        break;
                    case Value_NonInteger:
                        if(tree.IsAlwaysInteger()) return NoMatch;
                        break;
                }
                switch( ImmedConstraint_Sign(count & SignMask) )
                {
                    /*case SignMask: break;*/
                    case Sign_AnySign: break;
                    case Sign_Positive:
                        if(!tree.IsAlwaysSigned(true)) return NoMatch;
                        break;
                    case Sign_Negative:
                        if(!tree.IsAlwaysSigned(false)) return NoMatch;
                        break;
                    case Sign_NoIdea:
                        if(tree.IsAlwaysSigned(false)) return NoMatch;
                        if(tree.IsAlwaysSigned(true)) return NoMatch;
                        break;
                }
                switch( ImmedConstraint_Oneness(count & OnenessMask) )
                {
                    case OnenessMask: break;
                    case Oneness_Any: break;
                    case Oneness_One:    return NoMatch;
                    case Oneness_NotOne: return NoMatch;
                }

                return pack.flist[index].Match(tree, match, match_index);
            }
            default: // means groupfunction. No ImmedConstraint
            {
                if(!tree.IsImmed()) return NoMatch;
                double res = tree.GetImmed();
                if(transformation == Negate) res = -res;
                if(transformation == Invert) res = 1/res;
                double res2;
                if(!GetConst(match, res2)) return NoMatch;
                if(transf == Negate) res2 = -res2;
                if(transf == Invert) res2 = 1/res2;
                if(transf == NotThe) res2 = res2 != 0;
                /*std::cout << std::flush;
                fprintf(stderr, "Comparing %.20f and %.20f\n", res, res2);
                fflush(stderr);*/
                return FloatEqual(res, res2) ? FoundLastMatch : NoMatch;
            }
        }
        return NoMatch;
    }

    bool ParamSpec::GetConst(
        const MatchedParams::CodeTreeMatch& match,
        double& result) const
    {
        switch(OpcodeType(opcode))
        {
            case NumConstant:
                result = GetPackConst(index);
                break;
            case ImmedHolder:
            {
                std::map<unsigned, double>::const_iterator
                    i = match.ImmedMap.find(index);
                if(i == match.ImmedMap.end()) return false; // impossible
                result = i->second;
                //fprintf(stderr, "immedholder: %.20f\n", result);
                break;
            }
            case NamedHolder:
            {
                std::map<unsigned, MatchedParams::CodeTreeMatch::NamedItem>::const_iterator
                    i = match.NamedMap.find(index);
                if(i == match.NamedMap.end()) return false; // impossible
                result = (double) i->second.howmany;
                //fprintf(stderr, "namedholder: %.20f\n", result);
                break;
            }
            case RestHolder:
            {
                // Not enumerable
                return false;
            }
            case SubFunction:
            {
                // Not enumerable
                return false;
            }
            default:
            {
                switch(OPCODE(opcode))
                {
                    case cAdd:
                        result=0;
                        for(unsigned p=0; p<count; ++p)
                        {
                            double tmp;
                            if(!pack.plist[index+p].GetConst(match, tmp)) return false;
                            result += tmp;
                        }
                        break;
                    case cMul:
                        result=1;
                        for(unsigned p=0; p<count; ++p)
                        {
                            double tmp;
                            if(!pack.plist[index+p].GetConst(match, tmp)) return false;
                            result *= tmp;
                        }
                        break;
                    case cMin:
                        for(unsigned p=0; p<count; ++p)
                        {
                            double tmp;
                            if(!pack.plist[index+p].GetConst(match, tmp)) return false;
                            if(p == 0 || tmp < result) result = tmp;
                        }
                        break;
                    case cMax:
                        for(unsigned p=0; p<count; ++p)
                        {
                            double tmp;
                            if(!pack.plist[index+p].GetConst(match, tmp)) return false;
                            if(p == 0 || tmp > result) result = tmp;
                        }
                        break;
                    case cSin: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::sin(result); break;
                    case cCos: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::cos(result); break;
                    case cTan: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::tan(result); break;
                    case cAsin: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::asin(result); break;
                    case cAcos: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::acos(result); break;
                    case cAtan: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::atan(result); break;
                    case cSinh: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::sinh(result); break;
                    case cCosh: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::cosh(result); break;
                    case cTanh: if(!pack.plist[index].GetConst(match, result))return false;
                                 result = std::tanh(result); break;

                    case cAsinh: if(!pack.plist[index].GetConst(match, result))return false;
                                 result = fp_asinh(result); break;
                    case cAcosh: if(!pack.plist[index].GetConst(match, result))return false;
                                 result = fp_acosh(result); break;
                    case cAtanh: if(!pack.plist[index].GetConst(match, result))return false;
                                 result = fp_atanh(result); break;

                    case cCeil: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::ceil(result); break;
                    case cFloor: if(!pack.plist[index].GetConst(match, result))return false;
                                 result = std::floor(result); break;
                    case cLog: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::log(result); break;
                    case cExp: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::exp(result); break;
                    case cExp2: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::pow(2.0, result); break;
                    case cLog2: if(!pack.plist[index].GetConst(match, result))return false;
                                result = std::log(result) * CONSTANT_L2I;
                                //result = std::log2(result);
                                break;
                    case cLog10: if(!pack.plist[index].GetConst(match, result))return false;
                                 result = std::log10(result); break;
                    case cAbs: if(!pack.plist[index].GetConst(match, result))return false;
                               result = std::fabs(result); break;
                    case cPow:
                    {
                        if(!pack.plist[index+0].GetConst(match, result))return false;
                        double tmp;
                        if(!pack.plist[index+1].GetConst(match, tmp))return false;
                        result = std::pow(result, tmp);
                        //fprintf(stderr, "pow result: %.20f\n", result);
                        break;
                    }
                    case cMod:
                    {
                        if(!pack.plist[index+0].GetConst(match, result))return false;
                        double tmp;
                        if(!pack.plist[index+1].GetConst(match, tmp))return false;
                        result = std::fmod(result, tmp);
                        break;
                    }
                    default:
                        fprintf(stderr, "Unknown macro opcode: %s\n",
                            FP_GetOpcodeName(opcode).c_str());
                        return false;
                }
            }
        }
        if(transformation == Negate) result = -result;
        if(transformation == Invert) result = 1.0 / result;
        return true;
    }

    void MatchedParams::SynthesizeTree(
        FPoptimizer_CodeTree::CodeTree& tree,
        const MatchedParams& matcher,
        MatchedParams::CodeTreeMatch& match) const
    {
        for(unsigned a=0; a<count; ++a)
        {
            const ParamSpec& param = pack.plist[index+a];
            if(param.opcode == RestHolder)
            {
                // Add children directly to this tree
                param.SynthesizeTree(tree, matcher, match);
            }
            else
            {
                FPoptimizer_CodeTree::CodeTree* subtree = new FPoptimizer_CodeTree::CodeTree;
                param.SynthesizeTree(*subtree, matcher, match);
                subtree->ConstantFolding();
                subtree->Sort();
                subtree->Recalculate_Hash_NoRecursion(); // rehash this, but not the children, nor the parent
                FPoptimizer_CodeTree::CodeTree::Param p(subtree, param.sign) ;
                tree.AddParam(p);
            }
        }
    }

    void MatchedParams::ReplaceParams(
        FPoptimizer_CodeTree::CodeTree& tree,
        const MatchedParams& matcher,
        MatchedParams::CodeTreeMatch& match) const
    {
        // Replace the 0-level params indicated in "match" with the ones we have

        // First, construct the tree recursively using the "match" info
        SynthesizeTree(tree, matcher, match);

        // Remove the indicated params
        std::sort(match.param_numbers.begin(), match.param_numbers.end());
        for(size_t a=match.param_numbers.size(); a-->0; )
        {
            size_t num = match.param_numbers[a];
            tree.DelParam(num);
        }

        tree.ConstantFolding();

        tree.Sort();
        tree.Rehash(true); // rehash this and its parents, but not its children
    }

    void MatchedParams::ReplaceTree(
        FPoptimizer_CodeTree::CodeTree& tree,
        const MatchedParams& matcher,
        CodeTreeMatch& match) const
    {
        // Replace the entire tree with one indicated by our Params[0]
        // Note: The tree is still constructed using the holders indicated in "match".
        std::vector<FPoptimizer_CodeTree::CodeTree::Param> OldParams = tree.Params;
        tree.Params.clear();
        pack.plist[index].SynthesizeTree(tree, matcher, match);

        tree.ConstantFolding();

        tree.Sort();
        tree.Rehash(true);  // rehash this and its parents, but not its children
    }

    /* Synthesizes a new tree based on the given information
     * in ParamSpec. Assume the tree is empty, don't deallocate
     * anything. Don't touch Hash, Parent.
     */
    void ParamSpec::SynthesizeTree(
        FPoptimizer_CodeTree::CodeTree& tree,
        const MatchedParams& matcher,
        MatchedParams::CodeTreeMatch& match) const
    {
        switch(SpecialOpcode(opcode))
        {
            case RestHolder:
            {
                std::map<unsigned, std::vector<fphash_t> >
                    ::const_iterator i = match.RestMap.find(index);

                assert(i != match.RestMap.end());

                /*std::cout << std::flush;
                fprintf(stderr, "Restmap %u, sign %d, size is %u -- params %u\n",
                    (unsigned) i->first, sign, (unsigned) i->second.size(),
                    (unsigned) tree.Params.size());*/

                for(size_t a=0; a<i->second.size(); ++a)
                {
                    fphash_t hash = i->second[a];

                    std::map<fphash_t, FPoptimizer_CodeTree::CodeTreeP>
                        ::const_iterator j = match.trees.find(hash);

                    assert(j != match.trees.end());

                    FPoptimizer_CodeTree::CodeTree* subtree = j->second->Clone();
                    FPoptimizer_CodeTree::CodeTree::Param p(subtree, sign);
                    tree.AddParam(p);
                }
                /*fprintf(stderr, "- params size became %u\n", (unsigned)tree.Params.size());
                fflush(stderr);*/
                break;
            }
            case SubFunction:
            {
                const Function& fitem = pack.flist[index];
                tree.Opcode = fitem.opcode;
                const MatchedParams& mitem = pack.mlist[fitem.index];
                mitem.SynthesizeTree(tree, matcher, match);
                break;
            }
            case NamedHolder:
                if(!anyrepeat && minrepeat == 1)
                {
                    /* Literal parameter */
                    std::map<unsigned, MatchedParams::CodeTreeMatch::NamedItem>
                        ::iterator i = match.NamedMap.find(index);

                    assert(i != match.NamedMap.end());

                    fphash_t hash = i->second.hash;

                    std::map<fphash_t, FPoptimizer_CodeTree::CodeTreeP>
                        ::const_iterator j = match.trees.find(hash);

                    assert(j != match.trees.end());

                    tree.Opcode = j->second->Opcode;
                    switch(tree.Opcode)
                    {
                        case cImmed: tree.Value = j->second->Value; break;
                        case cVar:   tree.Var   = j->second->Var;  break;
                        case cFCall:
                        case cPCall: tree.Funcno = j->second->Funcno; break;
                    }

                    /* Note: SetParams() will Clone() all the given params.
                     *       This is considered appropriate, because the
                     *       same NamedHolder may be synthesized in multiple
                     *       trees.
                     *       Example of such rule:
                     *         asinh(x) -> log2(x + (x^2 + 1)^0.5) * CONSTANT_L2
                     *       We use n_synthesized here to limit the cloning only
                     *       to successive invokations of the same tree. The first
                     *       instance is simply assigned. This is safe, because the
                     *       tree from which it was brought, will not be used anymore.
                     */
                    tree.SetParams(j->second->Params, i->second.n_synthesized++ > 0);
                    break;
                }
                // passthru; x+ is synthesized as the number, not as the tree
            case NumConstant:
            case ImmedHolder:
            default:
                tree.Opcode = cImmed;
                GetConst(match, tree.Value); // note: return value is ignored
                // FIXME: Should we check ImmedConstraints here?
                break;
        }
    }

#ifdef DEBUG_SUBSTITUTIONS
    void DumpParam(const ParamSpec& p)
    {
        //std::cout << "/*p" << (&p-pack.plist) << "*/";

        if(p.sign) std::cout << '~';
        if(p.transformation == Negate) std::cout << '-';
        if(p.transformation == Invert) std::cout << '/';

        bool has_constraint = false;
        switch(SpecialOpcode(p.opcode))
        {
            case NumConstant: std::cout << GetPackConst(p.index); break;
            case ImmedHolder: has_constraint = true; std::cout << ImmedHolderNames[p.index]; break;
            case NamedHolder: has_constraint = true; std::cout << NamedHolderNames[p.index]; break;
            case RestHolder: std::cout << '<' << p.index << '>'; break;
            case SubFunction: DumpFunction(pack.flist[p.index]); break;
            default:
            {
                std::string opcode = FP_GetOpcodeName(p.opcode).substr(1);
                for(size_t a=0; a<opcode.size(); ++a) opcode[a] = std::toupper(opcode[a]);
                std::cout << opcode << '(';
                for(unsigned a=0; a<p.count; ++a)
                {
                    if(a > 0) std::cout << ' ';
                    DumpParam(pack.plist[p.index+a]);
                }
                std::cout << " )";
            }
        }
        if(p.anyrepeat && p.minrepeat==1) std::cout << '*';
        if(p.anyrepeat && p.minrepeat==2) std::cout << '+';
        if(has_constraint)
        {
            switch( ImmedConstraint_Value(p.count & ValueMask) )
            {
                case ValueMask: break;
                case Value_AnyNum: break;
                case Value_EvenInt:   std::cout << "@E"; break;
                case Value_OddInt:    std::cout << "@O"; break;
                case Value_IsInteger: std::cout << "@I"; break;
                case Value_NonInteger:std::cout << "@F"; break;
            }
            switch( ImmedConstraint_Sign(p.count & SignMask) )
            {
                case SignMask: break;
                case Sign_AnySign: break;
                case Sign_Positive:   std::cout << "@P"; break;
                case Sign_Negative:   std::cout << "@N"; break;
            }
            switch( ImmedConstraint_Oneness(p.count & OnenessMask) )
            {
                case OnenessMask: break;
                case Oneness_Any: break;
                case Oneness_One:     std::cout << "@1"; break;
                case Oneness_NotOne:  std::cout << "@M"; break;
            }
        }
    }

    void DumpParams(const MatchedParams& mitem)
    {
        //std::cout << "/*m" << (&mitem-pack.mlist) << "*/";

        if(mitem.type == PositionalParams) std::cout << '[';
        if(mitem.type == SelectedParams) std::cout << '{';

        for(unsigned a=0; a<mitem.count; ++a)
        {
            std::cout << ' ';
            DumpParam(pack.plist[mitem.index + a]);
        }

        switch(mitem.balance)
        {
            case BalanceMorePos: std::cout << " =+"; break;
            case BalanceMoreNeg: std::cout << " =-"; break;
            case BalanceEqual:   std::cout << " =="; break;
            case BalanceDontCare: break;
        }

        if(mitem.type == PositionalParams) std::cout << " ]";
        if(mitem.type == SelectedParams) std::cout << " }";
    }

    void DumpFunction(const Function& fitem)
    {
        //std::cout << "/*f" << (&fitem-pack.flist) << "*/";

        std::cout << '(' << FP_GetOpcodeName(fitem.opcode);
        DumpParams(pack.mlist[fitem.index]);
        std::cout << ')';
    }
    void DumpMatch(const Function& input,
                   const FPoptimizer_CodeTree::CodeTree& tree,
                   const MatchedParams& replacement,
                   const MatchedParams::CodeTreeMatch& matchrec,
                   bool DidMatch)
    {
        std::cout <<
            "Found " << (DidMatch ? "match" : "mismatch") << ":\n"
            "  Pattern    : ";
        DumpFunction(input);
        std::cout << "\n"
            "  Replacement: ";
        DumpParams(replacement);
        std::cout << "\n";

        std::cout <<
            "  Tree       : ";
        DumpTree(tree);
        std::cout << "\n";
        if(DidMatch) DumpHashes(tree);

        for(std::map<unsigned, MatchedParams::CodeTreeMatch::NamedItem>::const_iterator
            i = matchrec.NamedMap.begin(); i != matchrec.NamedMap.end(); ++i)
        {
            std::cout << "           " << NamedHolderNames[i->first] << " = ";
            DumpTree(*matchrec.trees.find(i->second.hash)->second);
            std::cout << " (" << i->second.howmany << " matches)\n";
        }

        for(std::map<unsigned, double>::const_iterator
            i = matchrec.ImmedMap.begin(); i != matchrec.ImmedMap.end(); ++i)
        {
            std::cout << "           " << ImmedHolderNames[i->first] << " = ";
            std::cout << i->second << std::endl;
        }

        for(std::map<unsigned, std::vector<fphash_t> >::const_iterator
            i = matchrec.RestMap.begin(); i != matchrec.RestMap.end(); ++i)
        {
            for(size_t a=0; a<i->second.size(); ++a)
            {
                fphash_t hash = i->second[a];
                std::cout << "         <" << i->first << "> = ";
                DumpTree(*matchrec.trees.find(hash)->second);
                std::cout << std::endl;
            }
            if(i->second.empty())
                std::cout << "         <" << i->first << "> = <empty>\n";
        }
        std::cout << std::flush;
    }
    void DumpHashes(const FPoptimizer_CodeTree::CodeTree& tree,
                    std::map<fphash_t, std::set<std::string> >& done)
    {
        for(size_t a=0; a<tree.Params.size(); ++a)
            DumpHashes(*tree.Params[a].param, done);

        std::stringstream buf;
        DumpTree(tree, buf);
        done[tree.Hash].insert(buf.str());
    }
    void DumpHashes(const FPoptimizer_CodeTree::CodeTree& tree)
    {
        std::map<fphash_t, std::set<std::string> > done;
        DumpHashes(tree, done);

        for(std::map<fphash_t, std::set<std::string> >::const_iterator
            i = done.begin();
            i != done.end();
            ++i)
        {
            const std::set<std::string>& flist = i->second;
            if(flist.size() != 1) std::cout << "ERROR - HASH COLLISION?\n";
            for(std::set<std::string>::const_iterator
                j = flist.begin();
                j != flist.end();
                ++j)
            {
                //std::cout << '[' << std::hex << i->first << ']' << std::dec;
                std::cout << ": " << *j << "\n";
            }
        }
    }
    void DumpTree(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o)
    {
        //o << "/*" << tree.Depth << "*/";
        const char* sep2 = "";
        //o << '[' << std::hex << tree.Hash << ']' << std::dec;
        switch(tree.Opcode)
        {
            case cImmed: o << tree.Value; return;
            case cVar:   o << "Var" << tree.Var; return;
            case cAdd: sep2 = " +"; break;
            case cMul: sep2 = " *"; break;
            case cAnd: sep2 = " &"; break;
            case cOr: sep2 = " |"; break;
            case cPow: sep2 = " ^"; break;
            default:
                o << FP_GetOpcodeName(tree.Opcode);
                if(tree.Opcode == cFCall || tree.Opcode == cPCall)
                    o << ':' << tree.Funcno;
        }
        o << '(';
        if(tree.Params.size() <= 1 && *sep2) o << (sep2+1) << ' ';
        for(size_t a=0; a<tree.Params.size(); ++a)
        {
            if(a > 0) o << ' ';
            if(tree.Params[a].sign) o << '~';

            DumpTree(*tree.Params[a].param, o);

            if(tree.Params[a].param->Parent != &tree)
            {
                o << "(?parent?)";
            }

            if(a+1 < tree.Params.size()) o << sep2;
        }
        o << ')';
    }
#endif
}

#endif
#include "io/fparser/fpconfig.h"
#include "io/fparser/fparser.h"
#include "io/fparser/fptypes.h"


using namespace FUNCTIONPARSERTYPES;

#ifdef FP_SUPPORT_OPTIMIZER
using namespace FPoptimizer_CodeTree;

void FunctionParser::Optimize()
{
    CopyOnWrite();

    //PrintByteCode(std::cout);

    FPoptimizer_CodeTree::CodeTreeP tree
        = CodeTree::GenerateFrom(data->ByteCode, data->Immed, *data);

    while(FPoptimizer_Grammar::pack.glist[0].ApplyTo(*tree)) // entry
        { //std::cout << "Rerunning 0\n";
        }

    while(FPoptimizer_Grammar::pack.glist[1].ApplyTo(*tree)) // intermediate
        { //std::cout << "Rerunning 1\n";
        }

    while(FPoptimizer_Grammar::pack.glist[2].ApplyTo(*tree)) // final1
        { //std::cout << "Rerunning 2\n";
        }

    while(FPoptimizer_Grammar::pack.glist[3].ApplyTo(*tree)) // final2
        { //std::cout << "Rerunning 3\n";
        }

    tree->Sort_Recursive();

    std::vector<unsigned> byteCode;
    std::vector<double> immed;
    size_t stacktop_max = 0;
    tree->SynthesizeByteCode(byteCode, immed, stacktop_max);

    /*std::cout << std::flush;
    std::cerr << std::flush;
    fprintf(stderr, "Estimated stacktop %u\n", (unsigned)stacktop_max);
    fflush(stderr);*/

    if(data->StackSize != stacktop_max)
    {
        data->StackSize = stacktop_max; // note: gcc warning is meaningful
        data->Stack.resize(stacktop_max);
    }

    data->ByteCode.swap(byteCode);
    data->Immed.swap(immed);

    //PrintByteCode(std::cout);
}

#endif
#include <cmath>
#include <list>
#include <cassert>

#include "io/fparser/fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

#ifndef FP_GENERATING_POWI_TABLE
static const unsigned MAX_POWI_BYTECODE_LENGTH = 15;
#else
static const unsigned MAX_POWI_BYTECODE_LENGTH = 999;
#endif
static const unsigned MAX_MULI_BYTECODE_LENGTH = 3;

//#define DEBUG_SUBSTITUTIONS

#ifdef DEBUG_SUBSTITUTIONS
namespace FPoptimizer_Grammar
{
    void DumpTree(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o = std::cout);
}
#endif

namespace
{
    using namespace FPoptimizer_CodeTree;

    bool AssembleSequence(
                  CodeTree& tree, long count,
                  const FPoptimizer_ByteCode::SequenceOpCode& sequencing,
                  FPoptimizer_ByteCode::ByteCodeSynth& synth,
                  size_t max_bytecode_grow_length);
}

namespace
{
    typedef
        std::map<fphash_t,  std::pair<size_t, CodeTreeP> >
        TreeCountType;
    typedef
        std::multimap<fphash_t, CodeTreeP>
        DoneTreesType;

    void FindTreeCounts(TreeCountType& TreeCounts, CodeTreeP tree)
    {
        TreeCountType::iterator i = TreeCounts.lower_bound(tree->Hash);
        if(i != TreeCounts.end()
        && tree->Hash == i->first
        && tree->IsIdenticalTo( * i->second.second ) )
            i->second.first += 1;
        else
            TreeCounts.insert(i, std::make_pair(tree->Hash, std::make_pair(size_t(1), tree)));

        for(size_t a=0; a<tree->Params.size(); ++a)
            FindTreeCounts(TreeCounts, tree->Params[a].param);
    }

    void RememberRecursivelyHashList(DoneTreesType& hashlist,
                                     const CodeTreeP& tree)
    {
        hashlist.insert( std::make_pair(tree->Hash, tree) );
        for(size_t a=0; a<tree->Params.size(); ++a)
            RememberRecursivelyHashList(hashlist, tree->Params[a].param);
    }
    void RecreateInversionsAndNegations(CodeTree& tree)
    {
        for(size_t a=0; a<tree.Params.size(); ++a)
            RecreateInversionsAndNegations(*tree.Params[a].param);

        bool changed = false;
        switch(tree.Opcode) // Recreate inversions and negations
        {
            case cMul:
            {
                for(size_t a=0; a<tree.Params.size(); ++a)
                    if(tree.Params[a].param->Opcode == cPow
                    && tree.Params[a].param->Params[1].param->IsImmed()
                    && tree.Params[a].param->Params[1].param->GetImmed() == -1)
                    {
                        tree.Params[a].param = tree.Params[a].param->Params[0].param;
                        tree.Params[a].param->Parent = &tree;
                        tree.Params[a].sign = !tree.Params[a].sign;
                        changed = true;
                    }
                break;
            }
            case cAdd:
            {
                for(size_t a=0; a<tree.Params.size(); ++a)
                    if(tree.Params[a].param->Opcode == cMul)
                    {
                        // if the mul group has a -1 constant...
                        bool subchanged = false;
                        CodeTree& mulgroup = *tree.Params[a].param;
                        for(size_t b=mulgroup.Params.size(); b-- > 0; )
                            if(mulgroup.Params[b].param->IsImmed()
                            && mulgroup.Params[b].param->GetImmed() == -1)
                            {
                                mulgroup.Params.erase(mulgroup.Params.begin()+b);
                                tree.Params[a].sign = !tree.Params[a].sign;
                                subchanged = true;
                            }
                        if(subchanged)
                        {
                            mulgroup.ConstantFolding();
                            mulgroup.Sort();
                            mulgroup.Recalculate_Hash_NoRecursion();
                            changed = true;
                        }
                    }
            }
        }
        if(changed)
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "BEGIN CONSTANTFOLDING: ";
            FPoptimizer_Grammar::DumpTree(tree);
            std::cout << "\n";
        #endif
            tree.ConstantFolding();
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "END CONSTANTFOLDING:   ";
            FPoptimizer_Grammar::DumpTree(tree);
            std::cout << "\n";
        #endif
            tree.Sort();
            tree.Rehash(true);
        }
    }
}

namespace FPoptimizer_CodeTree
{
    void CodeTree::SynthesizeByteCode(
        std::vector<unsigned>& ByteCode,
        std::vector<double>&   Immed,
        size_t& stacktop_max)
    {
        RecreateInversionsAndNegations(*this);

        FPoptimizer_ByteCode::ByteCodeSynth synth;

        /* Find common subtrees */
        TreeCountType TreeCounts;
        FindTreeCounts(TreeCounts, this);

        /* Synthesize some of the most common ones */
        DoneTreesType AlreadyDoneTrees;
    FindMore: ;
        size_t best_score = 0;
        TreeCountType::const_iterator synth_it;
        for(TreeCountType::const_iterator
            i = TreeCounts.begin();
            i != TreeCounts.end();
            ++i)
        {
            size_t score = i->second.first;
            // It must always occur at least twice
            if(score < 2) continue;
            // And it must not be a simple expression
            if(i->second.second->Depth < 2) CandSkip: continue;
            // And it must not yet have been synthesized
            DoneTreesType::const_iterator j = AlreadyDoneTrees.lower_bound(i->first);
            for(; j != AlreadyDoneTrees.end() && j->first == i->first; ++j)
            {
                if(j->second->IsIdenticalTo(*i->second.second))
                    goto CandSkip;
            }
            // Is a candidate.
            score *= i->second.second->Depth;
            if(score > best_score)
                { best_score = score; synth_it = i; }
        }
        if(best_score > 0)
        {
            /* Synthesize the selected tree */
            synth_it->second.second->SynthesizeByteCode(synth);
            /* Add the tree and all its children to the AlreadyDoneTrees list,
             * to prevent it from being re-synthesized
             */
            RememberRecursivelyHashList(AlreadyDoneTrees, synth_it->second.second);
            goto FindMore;
        }

        /* Then synthesize the actual expression */
        SynthesizeByteCode(synth);
      #ifndef FP_DISABLE_EVAL
        /* Ensure that the expression result is
         * the only thing that remains in the stack
         */
        /* Removed: Fparser does not seem to care! */
        /* But if cEval is supported, it still needs to be done. */
        if(synth.GetStackTop() > 1)
            synth.DoPopNMov(0, synth.GetStackTop()-1);
      #endif
        synth.Pull(ByteCode, Immed, stacktop_max);
    }

    void CodeTree::SynthesizeByteCode(FPoptimizer_ByteCode::ByteCodeSynth& synth)
    {
        // If the synth can already locate our operand in the stack,
        // never mind synthesizing it again, just dup it.
        /* FIXME: Possible hash collisions. */
        if(synth.FindAndDup(Hash))
        {
            return;
        }

        switch(Opcode)
        {
            case cVar:
                synth.PushVar(GetVar());
                break;
            case cImmed:
                synth.PushImmed(GetImmed());
                break;
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd:
            case cOr:
            {
                // Operand re-sorting:
                // If the first param has a sign, try to find a param
                // that does _not_ have a sign and put it first.
                // This can be done because params are commutative
                // when they are grouped with their signs.
                if(!Params.empty() && Params[0].sign)
                {
                    for(size_t a=1; a<Params.size(); ++a)
                        if(!Params[a].sign)
                        {
                            std::swap(Params[0], Params[a]);
                            break;
                        }
                }

                // Try to ensure that Immeds don't have a sign
                for(size_t a=0; a<Params.size(); ++a)
                {
                    CodeTreeP& param = Params[a].param;
                    if(Params[a].sign && param->IsImmed())
                        switch(Opcode)
                        {
                            case cAdd: param->NegateImmed(); Params[a].sign=false; break;
                            case cMul: if(param->GetImmed() == 0.0) break;
                                       param->InvertImmed(); Params[a].sign=false; break;
                            case cAnd:
                            case cOr:  param->NotTheImmed(); Params[a].sign=false; break;
                        }
                }

                if(Opcode == cMul) // Special treatment for cMul sequences
                {
                    // If the paramlist contains an Immed, and that Immed
                    // fits in a long-integer, try to synthesize it
                    // as add-sequences instead.
                    for(size_t a=0; a<Params.size(); ++a)
                    {
                        Param p = Params[a];
                        CodeTreeP& param = p.param;
                        if(!p.sign && param->IsLongIntegerImmed())
                        {
                            long value = param->GetLongIntegerImmed();
                            Params.erase(Params.begin()+a);

                            bool success = AssembleSequence(
                                *this, value, FPoptimizer_ByteCode::AddSequence,
                                synth,
                                MAX_MULI_BYTECODE_LENGTH);

                            // Readd the token so that we don't need
                            // to deal with allocationd/deallocation here.
                            Params.insert(Params.begin()+a, p);

                            if(success)
                            {
                                // this tree was treated just fine
                                synth.StackTopIs(Hash);
                                return;
                            }
                        }
                    }
                }

                int n_stacked = 0;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    CodeTreeP const & param = Params[a].param;
                    bool               sign = Params[a].sign;

                    param->SynthesizeByteCode(synth);
                    ++n_stacked;

                    if(sign) // Is the operand negated/inverted?
                    {
                        if(n_stacked == 1)
                        {
                            // Needs unary negation/invertion. Decide how to accomplish it.
                            switch(Opcode)
                            {
                                case cAdd:
                                    synth.AddOperation(cNeg, 1); // stack state: -1+1 = +0
                                    break;
                                case cMul:
                                    synth.AddOperation(cInv, 1); // stack state: -1+1 = +0
                                    break;
                                case cAnd:
                                case cOr:
                                    synth.AddOperation(cNot, 1); // stack state: -1+1 = +0
                                    break;
                            }
                            // Note: We could use RDiv or RSub when the first
                            // token is negated/inverted and the second is not, to
                            // avoid cNeg/cInv/cNot, but thanks to the operand
                            // re-sorting in the beginning of this code, this
                            // situation never arises.
                            // cNeg/cInv/cNot is only synthesized when the group
                            // consists entirely of negated/inverted items.
                        }
                        else
                        {
                            // Needs binary negation/invertion. Decide how to accomplish it.
                            switch(Opcode)
                            {
                                case cAdd:
                                    synth.AddOperation(cSub, 2); // stack state: -2+1 = -1
                                    break;
                                case cMul:
                                    synth.AddOperation(cDiv, 2); // stack state: -2+1 = -1
                                    break;
                                case cAnd:
                                case cOr:
                                    synth.AddOperation(cNot,   1);   // stack state: -1+1 = +0
                                    synth.AddOperation(Opcode, 2); // stack state: -2+1 = -1
                                    break;
                            }
                            n_stacked = n_stacked - 2 + 1;
                        }
                    }
                    else if(n_stacked > 1)
                    {
                        // Cumulate at the earliest opportunity.
                        synth.AddOperation(Opcode, 2); // stack state: -2+1 = -1
                        n_stacked = n_stacked - 2 + 1;
                    }
                }
                if(n_stacked == 0)
                {
                    // Uh, we got an empty cAdd/cMul/whatever...
                    // Synthesize a default value.
                    // This should never happen.
                    switch(Opcode)
                    {
                        case cAdd:
                        case cOr:
                            synth.PushImmed(0);
                            break;
                        case cMul:
                        case cAnd:
                            synth.PushImmed(1);
                            break;
                        case cMin:
                        case cMax:
                            //synth.PushImmed(NaN);
                            synth.PushImmed(0);
                            break;
                    }
                    ++n_stacked;
                }
                assert(n_stacked == 1);
                break;
            }
            case cPow:
            {
                const Param& p0 = Params[0];
                const Param& p1 = Params[1];

                if(p1.param->IsImmed() && p1.param->GetImmed() == 0.5)
                {
                    p0.param->SynthesizeByteCode(synth);
                    synth.AddOperation(cSqrt, 1);
                }
                else if(p1.param->IsImmed() && p1.param->GetImmed() == -0.5)
                {
                    p0.param->SynthesizeByteCode(synth);
                    synth.AddOperation(cRSqrt, 1);
                }
                /*
                else if(p0.param->IsImmed() && p0.param->GetImmed() == CONSTANT_E)
                {
                    p1.param->SynthesizeByteCode(synth);
                    synth.AddOperation(cExp, 1);
                }
                else if(p0.param->IsImmed() && p0.param->GetImmed() == CONSTANT_EI)
                {
                    p1.param->SynthesizeByteCode(synth);
                    synth.AddOperation(cNeg, 1);
                    synth.AddOperation(cExp, 1);
                }
                */
                else if(!p1.param->IsLongIntegerImmed()
                || !AssembleSequence( /* Optimize integer exponents */
                        *p0.param, p1.param->GetLongIntegerImmed(),
                        FPoptimizer_ByteCode::MulSequence,
                        synth,
                        MAX_POWI_BYTECODE_LENGTH)
                  )
                {
                    if(p0.param->IsImmed() && p0.param->GetImmed() > 0.0)
                    {
                        // Convert into cExp or Exp2.
                        //    x^y = exp(log(x) * y) =
                        //    Can only be done when x is positive, though.
                        double mulvalue = std::log( p0.param->GetImmed() );

                        if(p1.param->Opcode == cMul)
                        {
                            // Neat, we can delegate the multiplication to the child
                            p1.param->AddParam( Param(new CodeTree(mulvalue), false) );
                            p1.param->ConstantFolding();
                            p1.param->Sort();
                            p1.param->Recalculate_Hash_NoRecursion();
                            mulvalue = 1.0;
                        }

                        // If the exponent needs multiplication, multiply it
                        if(
                      #ifdef FP_EPSILON
                          fabs(mulvalue - (double)(long)mulvalue) <= FP_EPSILON
                      #else
                          mulvalue == (double)(long)mulvalue
                      #endif
                        && AssembleSequence(*p1.param, (long)mulvalue,
                                            FPoptimizer_ByteCode::AddSequence, synth,
                                            MAX_MULI_BYTECODE_LENGTH))
                        {
                            // Done with a dup/add sequence, cExp
                            synth.AddOperation(cExp, 1);
                        }
                        /* - disabled cExp2 optimizations for now, because it
                         *   turns out that glibc for at least x86_64 has a
                         *   particularly stupid exp2() implementation that
                         *   is _slower_ than exp() or even pow(2,x)
                         *
                        else if(
                          #ifndef FP_SUPPORT_EXP2
                           #ifdef FP_EPSILON
                            fabs(mulvalue - CONSTANT_L2) <= FP_EPSILON
                           #else
                            mulvalue == CONSTANT_L2
                           #endif
                          #else
                            true
                          #endif
                            )
                        {
                            // Do with cExp2; in all likelihood it's never slower than cExp.
                            mulvalue *= CONSTANT_L2I;
                            if(
                          #ifdef FP_EPSILON
                              fabs(mulvalue - (double)(long)mulvalue) <= FP_EPSILON
                          #else
                              mulvalue == (double)(long)mulvalue
                          #endif
                            && AssembleSequence(*p1.param, (long)mulvalue,
                                                FPoptimizer_ByteCode::AddSequence, synth,
                                                MAX_MULI_BYTECODE_LENGTH))
                            {
                                // Done with a dup/add sequence, cExp2
                                synth.AddOperation(cExp2, 1);
                            }
                            else
                            {
                                // Do with cMul and cExp2
                                p1.param->SynthesizeByteCode(synth);
                                synth.PushImmed(mulvalue);
                                synth.AddOperation(cMul, 2);
                                synth.AddOperation(cExp2, 1);
                            }
                        }*/
                        else
                        {
                            // Do with cMul and cExp
                            p1.param->SynthesizeByteCode(synth);
                            synth.PushImmed(mulvalue);
                            synth.AddOperation(cMul, 2);
                            synth.AddOperation(cExp, 1);
                        }
                    }
                    else
                    {
                        p0.param->SynthesizeByteCode(synth);
                        p1.param->SynthesizeByteCode(synth);
                        synth.AddOperation(Opcode, 2); // Create a vanilla cPow.
                    }
                }
                break;
            }
            case cIf:
            {
                size_t ofs;
                // If the parameter amount is != 3, we're screwed.
                Params[0].param->SynthesizeByteCode(synth); // expression
                synth.SynthIfStep1(ofs);
                Params[1].param->SynthesizeByteCode(synth); // true branch
                synth.SynthIfStep2(ofs);
                Params[2].param->SynthesizeByteCode(synth); // false branch
                synth.SynthIfStep3(ofs);
                break;
            }
            case cFCall:
            {
                // If the parameter count is invalid, we're screwed.
                for(size_t a=0; a<Params.size(); ++a)
                    Params[a].param->SynthesizeByteCode(synth);
                synth.AddOperation(Opcode, (unsigned) Params.size());
                synth.AddOperation(Funcno, 0, 0);
                break;
            }
            case cPCall:
            {
                // If the parameter count is invalid, we're screwed.
                for(size_t a=0; a<Params.size(); ++a)
                    Params[a].param->SynthesizeByteCode(synth);
                synth.AddOperation(Opcode, (unsigned) Params.size());
                synth.AddOperation(Funcno, 0, 0);
                break;
            }
            default:
            {
                // If the parameter count is invalid, we're screwed.
                for(size_t a=0; a<Params.size(); ++a)
                    Params[a].param->SynthesizeByteCode(synth);
                synth.AddOperation(Opcode, (unsigned) Params.size());
                break;
            }
        }
        synth.StackTopIs(Hash);
    }
}

namespace
{
    bool AssembleSequence(
        CodeTree& tree, long count,
        const FPoptimizer_ByteCode::SequenceOpCode& sequencing,
        FPoptimizer_ByteCode::ByteCodeSynth& synth,
        size_t max_bytecode_grow_length)
    {
        if(count != 0)
        {
            FPoptimizer_ByteCode::ByteCodeSynth backup = synth;

            tree.SynthesizeByteCode(synth);

            // Ignore the size generated by subtree
            size_t bytecodesize_backup = synth.GetByteCodeSize();

            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);

            size_t bytecode_grow_amount = synth.GetByteCodeSize() - bytecodesize_backup;
            if(bytecode_grow_amount > max_bytecode_grow_length)
            {
                synth = backup;
                return false;
            }
            return true;
        }
        else
        {
            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);
            return true;
        }
    }
}

#endif
#include <cmath>
#include <cassert>

#include "io/fparser/fptypes.h"
#include "io/fparser/fparser.h"


#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;


namespace FPoptimizer_CodeTree
{
    class CodeTreeParserData
    {
    private:
        std::vector<CodeTreeP> stack;
    public:
        CodeTreeParserData() : stack() { }

        void Eat(size_t nparams, OPCODE opcode)
        {
            CodeTreeP newnode = new CodeTree;
            newnode->Opcode = opcode;
            size_t stackhead = stack.size() - nparams;
            for(size_t a=0; a<nparams; ++a)
            {
                CodeTree::Param param;
                param.param = stack[stackhead + a];
                param.sign  = false;
                newnode->AddParam(param);
            }
            stack.resize(stackhead);
            stack.push_back(newnode);
        }

        void EatFunc(size_t params, OPCODE opcode, unsigned funcno)
        {
            Eat(params, opcode);
            stack.back()->Funcno = funcno;
        }

        void AddConst(double value)
        {
            CodeTreeP newnode = new CodeTree(value);
            stack.push_back(newnode);
        }

        void AddVar(unsigned varno)
        {
            CodeTreeP newnode = new CodeTree;
            newnode->Opcode = cVar;
            newnode->Var    = varno;
            stack.push_back(newnode);
        }

        void SwapLastTwoInStack()
        {
            std::swap(stack[stack.size()-1],
                      stack[stack.size()-2]);
        }

        void Dup()
        {
            Fetch(stack.size()-1);
        }

        void Fetch(size_t which)
        {
            stack.push_back(stack[which]->Clone());
        }

        void PopNMov(size_t target, size_t source)
        {
            stack[target] = stack[source];
            stack.resize(target+1);
        }

        CodeTreeP PullResult()
        {
            CodeTreeP result = stack.back();
            stack.resize(stack.size()-1);
            result->Rehash(false);
            result->Sort_Recursive();
            return result;
        }

        void CheckConst()
        {
            // Check if the last token on stack can be optimized with constant math
            CodeTreeP result = stack.back();
            result->ConstantFolding();
        }
    private:
        CodeTreeParserData(const CodeTreeParserData&);
        CodeTreeParserData& operator=(const CodeTreeParserData&);
    };

    CodeTreeP CodeTree::GenerateFrom(
        const std::vector<unsigned>& ByteCode,
        const std::vector<double>& Immed,
        const FunctionParser::Data& fpdata)
    {
        CodeTreeParserData data;
        std::vector<size_t> labels;

        for(size_t IP=0, DP=0; ; ++IP)
        {
            while(!labels.empty() && labels.back() == IP)
            {
                // The "else" of an "if" ends here
                data.Eat(3, cIf);
                labels.erase(labels.end()-1);
            }
            if(IP >= ByteCode.size()) break;

            unsigned opcode = ByteCode[IP];
            if(OPCODE(opcode) >= VarBegin)
            {
                data.AddVar(opcode);
            }
            else
            {
                switch( OPCODE(opcode) )
                {
                    // Specials
                    case cIf:
                        IP += 2;
                        continue;
                    case cJump:
                        labels.push_back(ByteCode[IP+1]+1);
                        IP += 2;
                        continue;
                    case cImmed:
                        data.AddConst(Immed[DP++]);
                        break;
                    case cDup:
                        data.Dup();
                        break;
                    case cNop:
                        break;
                    case cFCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        unsigned params = fpdata.FuncPtrs[funcno].params;
                        data.EatFunc(params, OPCODE(opcode), funcno);
                        break;
                    }
                    case cPCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        unsigned params = fpdata.FuncParsers[funcno].params;
                        data.EatFunc(params, OPCODE(opcode), funcno);
                        break;
                    }
                    // Unary operators requiring special attention
                    case cInv: // from fpoptimizer
                        data.AddConst(-1);
                        data.Eat(2, cPow); // 1/x is x^-1
                        break;
                    case cNeg:
                        data.AddConst(-1);
                        data.Eat(2, cMul); // -x is x*-1
                        break;
                    case cSqr: // from fpoptimizer
                        data.Dup();
                        data.Eat(2, cMul);
                        break;
                    // Unary functions requiring special attention
                    case cDeg:
                        data.AddConst(CONSTANT_DR);
                        data.Eat(2, cMul);
                        break;
                    case cRad:
                        data.AddConst(CONSTANT_RD);
                        data.Eat(2, cMul);
                        break;
                    case cExp:
                        data.AddConst(CONSTANT_E);
                        data.SwapLastTwoInStack();
                        data.Eat(2, cPow);
                        break;
                    case cExp2: // from fpoptimizer
                        data.AddConst(2.0);
                        data.SwapLastTwoInStack();
                        data.Eat(2, cPow);
                        break;
                    case cSqrt:
                        data.AddConst(0.5);
                        data.Eat(2, cPow);
                        break;
                    case cCot:
                        data.Eat(1, cTan);
                        data.AddConst(-1);
                        data.Eat(2, cPow);
                        break;
                    case cCsc:
                        data.Eat(1, cSin);
                        data.AddConst(-1);
                        data.Eat(2, cPow);
                        break;
                    case cSec:
                        data.Eat(1, cCos);
                        data.AddConst(-1);
                        data.Eat(2, cPow);
                        break;
                    case cLog10:
                        data.Eat(1, cLog);
                        data.AddConst(CONSTANT_L10I);
                        data.Eat(2, cMul);
                        break;
                    case cLog2:
                        data.Eat(1, cLog);
                        data.AddConst(CONSTANT_L2I);
                        data.Eat(2, cMul);
                        break;
                    //case cLog:
                    //    data.Eat(1, cLog2);
                    //    data.AddConst(CONSTANT_L2);
                    //    data.Eat(2, cMul);
                    //    break;
                    // Binary operators requiring special attention
                    case cSub:
                        data.AddConst(-1);
                        data.Eat(2, cMul); // -x is x*-1
                        data.Eat(2, cAdd); // Minus is negative adding
                        break;
                    case cRSub: // from fpoptimizer
                        data.SwapLastTwoInStack();
                        data.AddConst(-1);
                        data.Eat(2, cMul); // -x is x*-1
                        data.Eat(2, cAdd);
                        break;
                    case cDiv:
                        data.AddConst(-1);
                        data.Eat(2, cPow); // 1/x is x^-1
                        data.Eat(2, cMul); // Divide is inverse multiply
                        break;
                    case cRDiv: // from fpoptimizer
                        data.SwapLastTwoInStack();
                        data.AddConst(-1);
                        data.Eat(2, cPow); // 1/x is x^-1
                        data.Eat(2, cMul); // Divide is inverse multiply
                        break;
                    case cRSqrt: // from fpoptimizer
                        data.AddConst(-0.5);
                        data.Eat(2, cPow);
                        break;
                    // Binary operators not requiring special attention
                    case cAdd: case cMul:
                    case cMod: case cPow:
                    case cEqual: case cLess: case cGreater:
                    case cNEqual: case cLessOrEq: case cGreaterOrEq:
                    case cAnd: case cOr:
                        data.Eat(2, OPCODE(opcode));
                        break;
                    // Unary operators not requiring special attention
                    case cNot:
                    case cNotNot: // from fpoptimizer
                        data.Eat(1, OPCODE(opcode));
                        break;
                    // Special opcodes generated by fpoptimizer itself
                    case cFetch:
                        data.Fetch(ByteCode[++IP]);
                        break;
                    case cPopNMov:
                    {
                        unsigned stackOffs_target = ByteCode[++IP];
                        unsigned stackOffs_source = ByteCode[++IP];
                        data.PopNMov(stackOffs_target, stackOffs_source);
                        break;
                    }
                    // Note: cVar should never be encountered in bytecode.
                    // Other functions
#ifndef FP_DISABLE_EVAL
                    case cEval:
                    {
                        size_t paramcount = fpdata.variableRefs.size();
                        data.Eat(paramcount, OPCODE(opcode));
                        break;
                    }
#endif
                    default:
                        unsigned funcno = opcode-cAbs;
                        assert(funcno < FUNC_AMOUNT);
                        const FuncDefinition& func = Functions[funcno];
                        data.Eat(func.params, OPCODE(opcode));
                        break;
                }
            }
            data.CheckConst();
        }
        return data.PullResult();
    }
}

#endif

#endif
