# Patches

## GLFW color rendering on MacOs

Fix end of function `createNativeWindow` in `cocoa_window.mm` by adding

```
    [window->ns.object setColorSpace:[NSColorSpace sRGBColorSpace]];
```

Hit taken from 

```
    // https://github.com/google/filament/blob/main/libs/filamentapp/src/NativeWindowHelperCocoa.mm
    // function prepareNativeWindow
```
