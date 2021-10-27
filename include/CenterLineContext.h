#pragma once
namespace CenterLine{
    struct Context{
        bool prefer_ortho; // Orthogonal is preferred
        Context() : prefer_ortho(true) {}
    };
}
