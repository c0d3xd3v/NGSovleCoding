import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Item {
    SwipeView {
        id: swipeView
        anchors.fill: parent

        MainPage {
            id: mainPage
        }

        Loader {
            id: loader
            source: "RenderSettings.qml"
        }
    }
}
